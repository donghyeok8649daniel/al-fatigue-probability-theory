"""3D small-strain isotropic tetrahedral FEM with prescribed face tractions.

All six stress components and all assigned faces enter the mechanics. A pure
traction problem has six rigid motions: discrete orthogonality constraints select a
displacement representative, without introducing a physical clamp. Compatible
forces produce negligible gauge reactions. See SOLID_MECHANICS.md.
"""
from functools import lru_cache

import numpy as np
from scipy import sparse
from scipy.linalg import qr
from scipy.sparse.linalg import splu, cg, LinearOperator
from scipy.sparse.csgraph import connected_components

from .load_balance import (face_geometry, resultant, balanced_traction, balance_report,
                           LoadBalanceError, preflight_balance)
from .surface_setup import mesh_id
from .tensor_load import compile_tensor_matrix, evaluate_tensor

from .volume_mesher import MAX_SOLID_NODES, MAX_SOLID_CELLS, SolidMeshLimitError, tetrahedralize
ASSEMBLY_CHUNK_CELLS = 16_384
DIRECT_SOLVER_MAX_NODES = 80_000


def elastic_matrix(young_mpa, poisson):
    if not np.isfinite(young_mpa) or young_mpa <= 0 or not np.isfinite(poisson) or not -1 < poisson < .5:
        raise ValueError('solid.material')
    mu = young_mpa/(2*(1+poisson))
    lam = young_mpa*poisson/((1+poisson)*(1-2*poisson))
    c = np.zeros((6, 6))
    c[:3, :3] = lam
    c[np.arange(3), np.arange(3)] += 2*mu
    c[np.arange(3, 6), np.arange(3, 6)] = mu
    return c


def tensor_from_voigt(values):
    values = np.asarray(values)
    out = np.zeros((*values.shape[:-1], 3, 3))
    out[..., 0, 0], out[..., 1, 1], out[..., 2, 2] = np.moveaxis(values[..., :3], -1, 0)
    out[..., 1, 2] = out[..., 2, 1] = values[..., 3]
    out[..., 0, 2] = out[..., 2, 0] = values[..., 4]
    out[..., 0, 1] = out[..., 1, 0] = values[..., 5]
    return out


def von_mises(values):
    s = np.asarray(values)
    return np.sqrt(.5*((s[..., 0]-s[..., 1])**2+(s[..., 1]-s[..., 2])**2
                     +(s[..., 2]-s[..., 0])**2)+3*np.sum(s[..., 3:]**2, axis=-1))


def check_balance(mesh, traction, relative_tolerance=1e-9):
    report = balance_report(mesh, traction)
    net, relative = report['net'], report['relative']
    if relative > relative_tolerance:
        raise LoadBalanceError(net)
    return net, relative


def strain_matrix(gradients):
    """Expand the 12 gradient coefficients only for the current element chunk."""
    b = np.zeros((len(gradients), 6, 12))
    for j in range(4):
        gx, gy, gz = gradients[:, j].T
        b[:, 0, 3*j] = gx; b[:, 1, 3*j+1] = gy; b[:, 2, 3*j+2] = gz
        b[:, 3, 3*j+1] = gz; b[:, 3, 3*j+2] = gy
        b[:, 4, 3*j] = gz; b[:, 4, 3*j+2] = gx
        b[:, 5, 3*j] = gy; b[:, 5, 3*j+1] = gx
    return b


def element_dofs(cells):
    return (3*cells[..., None]+np.arange(3, dtype=np.int32)).reshape(-1, 12)


def element_stress(gradients, displacement, c):
    """Differentiate relative nodal motion before multiplying by stiffness.

    The four shape gradients sum to zero. Subtracting the first nodal motion
    enforces translation invariance before cancellation on very thin elements,
    and avoids expanding a 6x12 strain matrix for stress recovery.
    """
    difference = displacement[:, 1:]-displacement[:, :1]
    derivative = np.einsum('tni,tnj->tij', difference, gradients[:, 1:])
    strain = np.column_stack((derivative[:, 0, 0], derivative[:, 1, 1], derivative[:, 2, 2],
        derivative[:, 1, 2]+derivative[:, 2, 1], derivative[:, 0, 2]+derivative[:, 2, 0],
        derivative[:, 0, 1]+derivative[:, 1, 0]))
    return strain @ c.T


def assemble_stiffness(gradients, c, volumes, cells, size):
    """Bound dense work to one chunk and merge sparse chunks in a binary tree.

    A full 12x12 element array plus two 64-bit index arrays cost 3.5 GB at
    one million tetrahedra. Equal-sized sparse merges avoid repeatedly adding
    each small chunk to an ever-growing global matrix.
    """
    pending = []
    for start in range(0, len(volumes), ASSEMBLY_CHUNK_CELLS):
        end = start+ASSEMBLY_CHUNK_CELLS
        local = strain_matrix(gradients[start:end])
        ke = np.einsum('tki,kl,tlj->tij', local, c, local, optimize=True)
        ke *= volumes[start:end, None, None]
        ids = element_dofs(cells[start:end])
        rows = np.broadcast_to(ids[:, :, None], ke.shape).ravel()
        columns = np.broadcast_to(ids[:, None, :], ke.shape).ravel()
        part = sparse.coo_matrix((ke.ravel(), (rows, columns)), shape=(size, size)).tocsc()
        del ke, rows, columns
        level = 0
        while level < len(pending) and pending[level] is not None:
            part = pending[level]+part
            pending[level] = None
            level += 1
        if level == len(pending): pending.append(part)
        else: pending[level] = part
    result = None
    for part in pending:
        if part is not None: result = part if result is None else result+part
    return result


class SolidMechanics:
    """Reusable pure-traction stiffness; algebraic gauges add no physical support."""
    def __init__(self, mesh, young_mpa, poisson, target_mm, *, progress=None, stop_requested=None):
        if mesh is None:
            raise ValueError('solid.mesh')
        if not np.isfinite(target_mm) or target_mm <= 0:
            raise ValueError('solid.target')
        self.mesh, self.young_mpa, self.poisson = mesh, young_mpa, poisson
        self.progress = progress
        self.stop_requested = stop_requested or (lambda: False)
        self.c = elastic_matrix(young_mpa, poisson)
        # Use local coordinates throughout the stiffness assembly.
        self.origin = mesh.vertices.mean(axis=0)
        generated = tetrahedralize(mesh, target_mm, self.origin,
            max_nodes=MAX_SOLID_NODES, max_cells=MAX_SOLID_CELLS,
            progress=progress, stop_requested=self.stop_requested)
        for key, value in generated.items(): setattr(self, key, value)
        nodes, cells = self.nodes, self.cells
        if self.stop_requested(): raise InterruptedError('axial.cancelled')
        if progress: progress('solid.assembling', dict(nodes=len(nodes), cells=len(cells)))
        t = nodes[cells]
        jac = (t[:, 1:]-t[:, :1]).transpose(0, 2, 1)
        grad = np.linalg.inv(jac)
        self.gradients = np.concatenate((-grad.sum(axis=1)[:, None, :], grad), axis=1)
        n = 3*len(nodes)
        del grad, t, jac
        self.k = assemble_stiffness(self.gradients, self.c, self.volumes, self.cells, n)
        # Three star edges connect each tetrahedron just as all twelve directed
        # edges do; connected_components treats this graph as undirected.
        graph = sparse.coo_matrix((np.ones(len(cells)*3, dtype=np.uint8),
            (np.repeat(cells[:, 0], 3), cells[:, 1:].ravel())),
            shape=(len(nodes), len(nodes)))
        if connected_components(graph, directed=False, return_labels=False) != 1:
            raise ValueError('solid.disconnected')
        del graph
        rigid = np.zeros((len(nodes), 3, 6))
        rigid[:, :, :3] = np.eye(3)
        center = nodes.mean(axis=0)
        position = (nodes-center)/np.linalg.norm(np.ptp(nodes, axis=0))
        for j in range(3):
            rigid[:, :, j+3] = np.cross(np.eye(3)[j], position)
        self.rigid, _ = np.linalg.qr(rigid.reshape(n, 6))
        del rigid, position
        if self.stop_requested(): raise InterruptedError('axial.cancelled')
        self.iterative = len(nodes) > DIRECT_SOLVER_MAX_NODES
        if progress: progress('solid.iterative_setup' if self.iterative else 'solid.factorizing', {})
        self._prepare_solver()
        if self.stop_requested(): raise InterruptedError('axial.cancelled')

    @property
    def b(self):
        """Compatibility diagnostic; production expands B only in small chunks."""
        return strain_matrix(self.gradients)

    @property
    def dofs(self):
        return element_dofs(self.cells)

    def project(self, vector):
        return vector-self.rigid @ (self.rigid.T @ vector)

    def _prepare_solver(self):
        if self.iterative:
            try:
                import pyamg
            except ImportError as exc:
                raise ValueError('solid.amg_dependency') from exc
            # Rigid translations/rotations are exact null modes and also the
            # multigrid near-nullspace candidates. Project BOTH the operator
            # and preconditioner; no diagonal regularization or hidden support.
            matrix = (self.k/self.young_mpa).tobsr(blocksize=(3, 3))
            self.multilevel = pyamg.smoothed_aggregation_solver(matrix, B=self.rigid,
                symmetry='symmetric', smooth=('energy', {'krylov': 'cg'}),
                presmoother=('block_gauss_seidel', {'sweep': 'symmetric'}),
                postsmoother=('block_gauss_seidel', {'sweep': 'symmetric'}),
                max_coarse=96, coarse_solver='pinv')
            preconditioner = self.multilevel.aspreconditioner()
            self.operator = LinearOperator(matrix.shape,
                matvec=lambda v: self.project(matrix @ self.project(v)), dtype=float)
            self.preconditioner = LinearOperator(matrix.shape,
                matvec=lambda v: self.project(preconditioner @ self.project(v)), dtype=float)
            self.linear_solver = 'projected_CG_rigid_mode_smoothed_aggregation'
            return
        # Six independent rigid-mode rows select temporary gauge coordinates.
        # Project force and re-center displacement in solve(): this is exactly
        # the orthogonal Neumann gauge, including unbalanced auxiliary bases.
        # Dense global Lagrange-multiplier columns formerly caused large LU fill.
        _, _, pivots = qr(self.rigid.T, mode='economic', pivoting=True)
        self.gauge_dofs = pivots[:6]
        if np.linalg.matrix_rank(self.rigid[self.gauge_dofs]) != 6:
            raise ValueError('solid.singular')
        free = np.ones(self.k.shape[0], dtype=bool)
        free[self.gauge_dofs] = False
        self.free_dofs = np.flatnonzero(free)
        reduced = self.k[self.free_dofs][:, self.free_dofs]/self.young_mpa
        try:
            self.factor = splu(reduced, permc_spec='MMD_AT_PLUS_A',
                                options={'SymmetricMode': True})
        except RuntimeError as exc:
            raise ValueError('solid.singular') from exc
        self.linear_solver = 'sparse_LU_equivalent_orthogonal_Neumann_gauge'

    def nodal_force(self, traction):
        forces = np.zeros((len(self.nodes), 3))
        child_force = np.asarray(traction)[self.parents]*self.boundary_areas[:, None]/3
        for i in range(3):
            np.add.at(forces, self.boundary_faces[:, i], child_force)
        return forces.ravel()

    def solve(self, traction, *, require_balance=True):
        if require_balance:
            check_balance(self.mesh, traction)
        force = self.nodal_force(traction)
        gauge_reaction = self.rigid @ (self.rigid.T @ force)
        compatible = force-gauge_reaction
        iterations = 0
        if self.stop_requested(): raise InterruptedError('axial.cancelled')
        if self.iterative:
            def iteration(_u):
                nonlocal iterations
                iterations += 1
                if self.stop_requested(): raise InterruptedError('axial.cancelled')
                if self.progress and iterations % 25 == 0:
                    self.progress('solid.iterating', dict(iterations=iterations))
            u, info = cg(self.operator, compatible/self.young_mpa, M=self.preconditioner,
                         rtol=1e-11, atol=0., maxiter=2000, callback=iteration)
            if info != 0: raise ValueError('solid.iterative_failed')
        else:
            u = np.zeros_like(force)
            u[self.free_dofs] = self.factor.solve(compatible[self.free_dofs]/self.young_mpa)
        u = self.project(u)
        residual = self.k @ u-force
        stress = np.empty((len(self.cells), 6))
        for start in range(0, len(self.cells), ASSEMBLY_CHUNK_CELLS):
            end = start+ASSEMBLY_CHUNK_CELLS
            stress[start:end] = element_stress(self.gradients[start:end],
                u.reshape(-1, 3)[self.cells[start:end]], self.c)
        face_stress = np.zeros((len(self.mesh.faces), 6))
        np.add.at(face_stress, self.parents, stress[self.boundary_cells]*self.boundary_areas[:, None])
        face_stress /= self.mesh.areas[:, None]
        relative = np.linalg.norm(residual)/np.linalg.norm(force) if np.any(force) else np.linalg.norm(residual)
        if require_balance and (not np.isfinite(relative) or relative > 1e-7):
            raise ValueError('solid.residual')
        return dict(displacement_mm=u.reshape(-1, 3), stress_voigt_mpa=stress,
            face_stress_voigt_mpa=face_stress, nodal_force_n=force.reshape(-1, 3),
            residual_n=residual.reshape(-1, 3), relative_residual=float(relative),
            linear_iterations=iterations,
            gauge_reaction_n=gauge_reaction.reshape(-1, 3),
            strain_energy_n_mm=float(.5*u @ (self.k @ u)))


class FaceStressHistory:
    """Exact linear mechanical superposition with actual time-expression evaluation.

    Individual basis loads need not balance; their gauge forces are auxiliary.
    The combined load is checked before every requested physical evaluation.
    An already approved boundary correction is applied to every basis, so its
    full time dependence is retained. No unapproved correction is inserted.
    """
    def __init__(self, system, loads, correction=None):
        self.system, self.loads = system, tuple(loads)
        self.compiled = [compile_tensor_matrix(load.tensor_expression) for load in loads]
        _, normals = face_geometry(system.mesh)
        size = 6*len(loads)
        self.traction_basis = np.empty((size, len(system.mesh.faces), 3))
        self.face_basis = np.empty((size, len(system.mesh.faces), 6))
        self.cell_basis = np.empty((size, len(system.cells), 6))
        self.displacement_basis = np.empty((size, len(system.nodes), 3))
        self.residual_basis = np.empty_like(self.displacement_basis)
        self.force_basis = np.empty_like(self.displacement_basis)
        self.linear_iterations = np.zeros(size, dtype=int)
        for i, load in enumerate(loads):
            ids = np.asarray(load.face_indices, dtype=int)
            if not len(ids) or ids.min() < 0 or ids.max() >= len(system.mesh.faces):
                raise ValueError('solid.loads')
            for j in range(6):
                if system.progress: system.progress('solid.basis_progress', dict(index=6*i+j+1, total=size))
                unit = tensor_from_voigt(np.eye(6)[j])
                traction = np.zeros((len(system.mesh.faces), 3))
                traction[ids] = normals[ids] @ unit.T
                if correction is not None:
                    traction = balanced_traction(system.mesh, traction, correction)
                result = system.solve(traction, require_balance=False)
                index = 6*i+j
                self.traction_basis[index] = traction
                self.face_basis[index] = result['face_stress_voigt_mpa']
                self.cell_basis[index] = result['stress_voigt_mpa']
                self.displacement_basis[index] = result['displacement_mm']
                self.residual_basis[index] = result['residual_n']
                self.force_basis[index] = result['nodal_force_n']
                self.linear_iterations[index] = result['linear_iterations']
                del result
        self.balance_centers, _ = face_geometry(system.mesh)
        self.balance_centers -= system.mesh.vertices.mean(axis=0)
        self.balance_length = np.linalg.norm(np.ptp(system.mesh.vertices, axis=0))
        self.surface_areas = system.mesh.areas
        flat_force = self.force_basis.reshape(size, 3*len(system.nodes))
        self.force_gram = flat_force @ flat_force.T
        self.force_basis_norm = np.linalg.norm(flat_force, axis=1)
        flat_residual = self.residual_basis.reshape(size, 3*len(system.nodes))
        self.rigid_residual = flat_residual @ system.rigid
        remainder = flat_residual-self.rigid_residual @ system.rigid.T
        self.remainder_gram = remainder @ remainder.T
        self._cached_coefficients = lru_cache(maxsize=128)(self._coefficients)

    def coefficients(self, time):
        return self._cached_coefficients(float(time))

    def _coefficients(self, time):
        coefficients = []
        for load, compiled in zip(self.loads, self.compiled):
            tensor = evaluate_tensor(compiled, t=time, frequency=load.frequency,
                normal_mean=load.normal_mean_mpa, normal_amplitude=load.normal_amplitude_mpa,
                shear_mean=load.shear_mean_mpa, shear_amplitude=load.shear_amplitude_mpa)
            if not np.allclose(tensor, tensor.T, rtol=1e-12, atol=1e-12):
                raise ValueError('solid.traction')
            coefficients.extend(tensor[[0, 1, 2, 1, 0, 0], [0, 1, 2, 2, 2, 1]])
        coefficients = np.array(coefficients)
        traction = np.einsum('b,bfi->fi', coefficients, self.traction_basis)
        # Static geometry and small Gram matrices avoid meshing/volume audits
        # and full nodal vector construction at every adaptive PDE time step.
        forces = traction*self.surface_areas[:, None]
        net = np.r_[forces.sum(axis=0), np.cross(self.balance_centers, forces).sum(axis=0)]
        scale = float(np.linalg.norm(forces, axis=1).sum())
        relative = float(np.max(np.abs(np.r_[net[:3], net[3:]/self.balance_length]))/scale) if scale else 0.
        if relative > 1e-9:
            raise LoadBalanceError(net, time)
        force_square = float(coefficients @ self.force_gram @ coefficients)
        if force_square <= 1e-12*(np.abs(coefficients) @ self.force_basis_norm)**2:
            force_square = float(np.sum(np.einsum('b,bni->ni', coefficients, self.force_basis)**2))
        residual_square = float(coefficients @ self.remainder_gram @ coefficients)
        residual_square += float(np.sum((coefficients @ self.rigid_residual)**2))
        if residual_square < 0:
            residual_square = float(np.sum(np.einsum('b,bni->ni', coefficients, self.residual_basis)**2))
        linear_relative = np.sqrt(residual_square/force_square) if force_square else np.sqrt(residual_square)
        if not np.isfinite(linear_relative) or linear_relative > 1e-7:
            raise ValueError('solid.residual')
        return coefficients, net, relative, float(linear_relative)

    def face_stress(self, time):
        coefficients, _, _, _ = self.coefficients(time)
        return np.einsum('b,bfc->fc', coefficients, self.face_basis)


def select_stress_histories(histories, count):
    """Deterministic farthest-point sampling in signed stress-history space.

    Actual PDEs run at the selected cells. Other cells display the nearest
    representative's probability. The sampled stress mismatch is retained;
    this finite spatial sampling is explicitly uncertified, not interpolation
    of P by triangle area or a manufactured stress concentration.
    """
    histories = np.asarray(histories)
    if int(count) != count or not 2 <= count <= 64:
        raise ValueError('axial.count_error')
    chosen = [int(np.argmax(np.max(np.abs(histories), axis=0)))]
    difference = np.max(np.abs(histories-histories[:, chosen[0], None]), axis=0)
    membership = np.zeros(histories.shape[1], dtype=int)
    scale = max(float(np.max(np.abs(histories))), np.finfo(float).tiny)
    while len(chosen) < min(count, histories.shape[1]) and difference.max() > 1e-10*scale:
        candidate = int(np.argmax(difference))
        distance = np.max(np.abs(histories-histories[:, candidate, None]), axis=0)
        closer = distance < difference
        membership[closer] = len(chosen)
        difference[closer] = distance[closer]
        chosen.append(candidate)
    return np.array(chosen), membership, difference


def run_solid_probability(config, mesh, loads, correction, *, poisson, target_mm,
                          direction, sample_count=16, stop_requested=None, progress=None,
                          record_callback=None):
    """3D tensor mechanics plus explicitly axial-projected local reference PDEs.

    The full tensor is retained and available in the map. The canonical local
    model receives e.T@sigma@e; this is a declared projection, not multiaxial
    fatigue calibration or a replacement of the underlying energy.
    """
    from .solver_adapter import run_ui_analysis
    from solver_v1.kinetic_calibration_workflow import build_time_basis_model
    from .surface_setup import encode
    config.require_material_backend()
    config.validate()
    stopped = stop_requested or (lambda: False)
    axis = np.asarray(direction, dtype=float)
    if axis.shape != (3,) or not np.isfinite(axis).all() or not np.linalg.norm(axis):
        raise ValueError('error.direction')
    axis = axis/np.linalg.norm(axis)
    if config.time_basis != 'model':
        raise ValueError('solid.time')
    if mesh is None:
        raise ValueError('solid.mesh')
    setup_json = encode(mesh, loads, correction)
    if stopped():
        raise InterruptedError('axial.cancelled')
    if progress: progress('solid.checking_loads', {})
    preflight_balance(mesh, loads, correction, config.cycles*config.model_period)
    if progress:
        progress('solid.meshing', {})
    system = SolidMechanics(mesh, config.young_mpa, poisson, target_mm,
                            progress=progress, stop_requested=stopped)
    history = FaceStressHistory(system, loads, correction)
    x, y, z = axis
    projection = np.array([x*x, y*y, z*z, 2*y*z, 2*x*z, 2*x*y])
    projected_basis = history.cell_basis @ projection
    end = config.cycles*config.model_period
    interval = config.model_period/config.steps_per_cycle
    times = np.arange(int(np.floor(config.cycles*config.steps_per_cycle))+1)*interval
    if end-times[-1] > 1e-12*config.model_period:
        times = np.append(times, end)
    else:
        times[-1] = end
    coefficients, nets, balances, residuals = [], [], [], []
    for time in times:
        if stopped():
            raise InterruptedError('axial.cancelled')
        coef, net, balance, residual = history.coefficients(float(time))
        coefficients.append(coef); nets.append(net); balances.append(balance); residuals.append(residual)
    coefficients = np.array(coefficients)
    # Bounded sampling matrix; all PDE-time evaluations still use the original
    # time expressions and check force/torque compatibility at the evaluated time.
    selection_times = np.linspace(0., end, min(max(len(times), 17), 129))
    selection_coefficients = np.array([history.coefficients(float(t))[0] for t in selection_times])
    samples, mapping, mismatch = select_stress_histories(selection_coefficients @ projected_basis, sample_count)
    raw_histories, probabilities, floors = [], [], []
    probability_model = build_time_basis_model(config.energy_model, time_basis=config.time_basis,
                                              calibration=config.time_calibration)
    for i, cell in enumerate(samples):
        if stopped():
            raise InterruptedError('axial.cancelled')
        if progress:
            progress('solid.progress', dict(index=i+1, total=len(samples)))
        def stress(time, cell=cell):
            return float(history.coefficients(float(time))[0] @ projected_basis[:, cell])
        # Main-window histories use the first representative PDE, exactly as the
        # returned reference_result does. Do not concatenate other cells' clocks
        # or present a partly computed spatial field as a finished specimen.
        local = run_ui_analysis(config, axial_stress_function=stress, stop_requested=stopped,
                                _prepared_model=probability_model,
                                record_callback=record_callback if i == 0 else None)
        if i == 0:
            reference_result = local
        if stopped():
            raise InterruptedError('axial.cancelled')
        actual_time = np.asarray(local['model_time'])
        extinct = float(np.asarray(local['intact_probability_mass'])[-1]) == 0.
        if abs(actual_time[-1]-end) > 1e-9*config.model_period and not extinct:
            raise ValueError('axial.incomplete')
        raw = {key: np.asarray(local[key]) for key in ('model_time', 'local_initiation_probability',
                'local_rare_event_floor', 'mass_balance_residual', 'applied_stress_mpa')}
        raw['cell_index'] = int(cell)
        raw['constant_extension_after_extinction'] = extinct and actual_time[-1] < end
        raw_histories.append(raw)
        probabilities.append(np.interp(times, actual_time, raw['local_initiation_probability']))
        floors.append(np.full_like(times, np.max(raw['local_rare_event_floor'])))
    p = np.array(probabilities).T
    return reference_result, dict(schema='aft.solid-neumann/1', mesh_sha256=mesh_id(mesh),
        mechanics='3D_linear_isotropic_tetrahedral_FEM_pure_traction',
        linear_solver=system.linear_solver,
        meshing_info=system.meshing_info,
        linear_basis_iterations=history.linear_iterations,
        young_mpa=config.young_mpa, poisson_ratio=poisson, target_mm=target_mm,
        nodes_mm=system.nodes+system.origin, tetrahedra=system.cells,
        cell_volumes_mm3=system.volumes, face_stress_basis_mpa=history.face_basis,
        cell_stress_basis_mpa=history.cell_basis, displacement_basis_mm=history.displacement_basis,
        load_coefficients=coefficients, model_time=times,
        force_torque_balance=np.array(nets), relative_balance=np.array(balances),
        relative_linear_residual=np.array(residuals),
        gauge='six_orthogonal_rigid_motion_constraints_no_physical_clamp',
        setup_json=setup_json,
        probability_projection='signed_axial_e_T_sigma_e_into_canonical_reduced_reference',
        probability_axis=axis, sample_cells=samples, cell_samples=mapping,
        sample_stress_history_mpa=coefficients @ projected_basis[:, samples],
        sampling_check_times=selection_times, sampling_stress_mismatch_mpa=mismatch,
        spatial_sampling='nearest_representative_in_sampled_signed_stress_history',
        local_initiation_probability=p, local_survival_probability=1-p,
        local_rare_event_floor=np.array(floors).T, raw_histories=raw_histories,
        probability_resolution_certified=False, spatial_resolution_certified=False,
        material_calibrated=False, energy_model=config.energy_model,
        energy_model_metadata={key: value for key, value in reference_result.items()
                               if key.startswith('energy_model_')},
        analysis_quality=config.analysis_quality, integration_method=config.integration_method,
        initialization=config.initialization,
        grid_shape=(config.grid_n_a, config.grid_n_s), time_basis=config.time_basis,
        units=dict(position='mm', stress='MPa', force='N', torque='N mm', probability='1'))


def solid_field(result, field, index):
    """Render stored data only; never rerun mechanics, time expressions or PDE."""
    if field in ('initiation', 'survival'):
        key = 'local_initiation_probability' if field == 'initiation' else 'local_survival_probability'
        return np.asarray(result[key])[index, np.asarray(result['cell_samples'])]
    coefficient, basis = result['load_coefficients'][index], result['cell_stress_basis_mpa']
    if field == 'stress':
        x, y, z = result['probability_axis']
        projection = np.array([x*x, y*y, z*z, 2*y*z, 2*x*z, 2*x*y])
        # Read only the required components; do not expand a 3x3 tensor per cell.
        value = np.zeros(basis.shape[1])
        for component in np.flatnonzero(projection):
            value += projection[component]*np.einsum('b,bt->t', coefficient, basis[..., component])
        return value
    if field == 'mises':
        return von_mises(np.einsum('b,btc->tc', coefficient, basis))
    component = {'xx': 0, 'yy': 1, 'zz': 2, 'yz': 3, 'xz': 4, 'xy': 5}[field]
    return np.einsum('b,bt->t', coefficient, basis[..., component])


def validate_solid_result(result, mesh=None):
    if not isinstance(result, dict) or result.get('schema') != 'aft.solid-neumann/1':
        raise ValueError('axial.incomplete')
    times, basis, coef = (np.asarray(result[k]) for k in ('model_time', 'cell_stress_basis_mpa', 'load_coefficients'))
    p, mapping = np.asarray(result['local_initiation_probability']), np.asarray(result['cell_samples'])
    if (times.ndim != 1 or not len(times) or not np.isfinite(times).all() or np.any(np.diff(times) <= 0)
            or basis.ndim != 3 or basis.shape[-1] != 6 or not np.isfinite(basis).all()
            or coef.shape != (len(times), len(basis)) or not np.isfinite(coef).all()
            or p.ndim != 2 or p.shape[0] != len(times) or not p.shape[1] or not np.isfinite(p).all()
            or np.any((p < 0) | (p > 1)) or np.any(np.diff(p, axis=0) < 0)
            or mapping.shape != (basis.shape[1],) or mapping.dtype.kind not in 'iu'
            or not len(mapping) or mapping.min() < 0 or mapping.max() >= p.shape[1]):
        raise ValueError('axial.incomplete')
    nodes, cells, volumes = (np.asarray(result[k]) for k in
                            ('nodes_mm', 'tetrahedra', 'cell_volumes_mm3'))
    axis, samples = np.asarray(result['probability_axis']), np.asarray(result['sample_cells'])
    survival, floor = (np.asarray(result[k]) for k in
                       ('local_survival_probability', 'local_rare_event_floor'))
    if (nodes.ndim != 2 or nodes.shape[1] != 3 or len(nodes) < 4 or not np.isfinite(nodes).all()
            or cells.shape != (len(mapping), 4) or cells.dtype.kind not in 'iu'
            or cells.min() < 0 or cells.max() >= len(nodes)
            or volumes.shape != (len(cells),) or not np.isfinite(volumes).all() or np.any(volumes <= 0)
            or axis.shape != (3,) or not np.isfinite(axis).all() or not np.isclose(np.linalg.norm(axis), 1.)
            or samples.shape != (p.shape[1],) or samples.dtype.kind not in 'iu'
            or samples.min() < 0 or samples.max() >= len(cells) or len(np.unique(samples)) != len(samples)
            or survival.shape != p.shape or not np.allclose(survival, 1-p, rtol=0, atol=1e-14)
            or floor.shape != p.shape or not np.isfinite(floor).all() or np.any(floor < 0)):
        raise ValueError('axial.incomplete')
    tets = nodes[cells]
    actual_volumes = np.abs(np.linalg.det((tets[:, 1:]-tets[:, :1]).transpose(0, 2, 1)))/6
    if not np.allclose(actual_volumes, volumes, rtol=1e-8, atol=0):
        raise ValueError('axial.incomplete')
    for key, shape in (
        ('force_torque_balance', (len(times), 6)),
        ('relative_balance', (len(times),)),
        ('relative_linear_residual', (len(times),)),
        ('sampling_stress_mismatch_mpa', (len(cells),)),
        ('sample_stress_history_mpa', p.shape),
        ('displacement_basis_mm', (len(basis), len(nodes), 3)),
    ):
        values = np.asarray(result[key])
        if values.shape != shape or not np.isfinite(values).all():
            raise ValueError('axial.incomplete')
    raw = result['raw_histories']
    if len(raw) != len(samples):
        raise ValueError('axial.incomplete')
    for sample, history in zip(samples, raw):
        raw_times = np.asarray(history['model_time'])
        if (history['cell_index'] != sample or raw_times.ndim != 1 or not len(raw_times)
                or not np.isfinite(raw_times).all() or np.any(np.diff(raw_times) <= 0)):
            raise ValueError('axial.incomplete')
        for key in ('local_initiation_probability', 'local_rare_event_floor',
                    'mass_balance_residual', 'applied_stress_mpa'):
            values = np.asarray(history[key])
            allowed_shape = values.shape == raw_times.shape or (key == 'local_rare_event_floor' and values.ndim == 0)
            if not allowed_shape or not np.isfinite(values).all():
                raise ValueError('axial.incomplete')
    if mesh is not None and result['mesh_sha256'] != mesh_id(mesh):
        raise ValueError('axial.stale')


def visible_cell_surface(result, fraction=1.):
    """Boundary of the displayed tetrahedron subset, plus each face's cell ID.

    This is visualization clipping by cell center, not a new mechanics solve.
    Internal faces cancel by vertex IDs; cut faces reveal the interior values.
    """
    nodes, cells = np.asarray(result['nodes_mm']), np.asarray(result['tetrahedra'])
    if not np.isfinite(fraction) or not 0 <= fraction <= 1:
        raise ValueError('invalid slice fraction')
    centers = nodes[cells].mean(axis=1)
    depth = centers @ np.asarray(result['probability_axis'])
    threshold = depth.min()+fraction*np.ptp(depth)
    ids = np.flatnonzero(depth <= threshold)
    selected = cells[ids]
    faces = np.concatenate([selected[:, [0, 2, 1]], selected[:, [0, 1, 3]],
                            selected[:, [0, 3, 2]], selected[:, [1, 2, 3]]])
    _, first, counts = np.unique(np.sort(faces, axis=1), axis=0, return_index=True, return_counts=True)
    boundary = first[counts == 1]
    return faces[boundary], np.tile(ids, 4)[boundary]
