"""Spatial original-SW crack reference, with no prescribed failure strength.

The energy is the SAME sum of atomic site energies as the rigid Si control.
Atoms can move independently; an opened initial crack is a configuration, not
a list of disabled bonds. Neighbors are rebuilt from geometry, including every
affected exterior site's angular energy. Units: Angstrom, eV, eV/Angstrom.

This is a 0 K constrained-energy research tool. Its relaxed energy is not a
finite-temperature PMF, and no physical mobility, time or wafer strength is
supplied. No production energy registry or probability operator is modified.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product

import numpy as np
from scipy.optimize import minimize
from scipy.spatial import cKDTree

from .silicon_environment_research import _parameters, environment_jet


@dataclass
class AtomicEvaluation:
    energy: float
    gradient: np.ndarray
    site_energy: np.ndarray


class SpatialSW:
    """Exact finite SW reference; optional periodicity along Cartesian z only.

    The Verlet list is valid only while EVERY atom is within skin/2 of its
    construction position. Periodic images are explicit, so thin front cells
    do not silently lose multiple images of a neighbor. The skin changes the
    candidate list only, never the physical SW cutoff.
    """
    def __init__(self, parameters, *, front_period=None, skin=.6):
        self.p = _parameters(parameters)
        if front_period is not None and (not np.isfinite(front_period) or front_period <= 0):
            raise ValueError('positive finite front period required')
        if not np.isfinite(skin) or skin < 0:
            raise ValueError('nonnegative finite neighbor skin required')
        self.period, self.skin = front_period, float(skin)
        self.cutoff = self.p['a']*self.p['sigma']
        self._reference = None
        self.rebuilds = 0

    def _edges(self, positions):
        r = np.asarray(positions, float)
        if r.ndim != 2 or r.shape[1] != 3 or not len(r) or not np.all(np.isfinite(r)):
            raise ValueError('nonempty finite (atoms,3) positions required')
        if (self._reference is None or r.shape != self._reference.shape
                or np.max(np.linalg.norm(r-self._reference, axis=1)) >= self.skin/2):
            wrapped = r.copy()
            wraps = np.zeros(len(r), int)
            if self.period is not None:
                wraps = np.floor(r[:, 2]/self.period).astype(int)
                wrapped[:, 2] -= wraps*self.period
                shell = int(np.ceil((self.cutoff+self.skin)/self.period))
                images = np.arange(-shell, shell+1)
            else:
                images = np.array([0])
            replicas = np.tile(wrapped, (len(images), 1))
            if self.period is not None:
                replicas[:, 2] += np.repeat(images*self.period, len(r))
            lists = cKDTree(replicas).query_ball_point(wrapped, self.cutoff+self.skin)
            i = np.repeat(np.arange(len(r)), [len(x) for x in lists])
            ids = np.concatenate(lists).astype(int)
            j, shifts = ids % len(r), images[ids//len(r)]
            keep = (i != j) | (shifts != 0)
            i, j, shifts = i[keep], j[keep], shifts[keep]
            self._i, self._j = i, j
            self._shifts = shifts+wraps[i]-wraps[j]
            self._reference = r.copy()
            self.rebuilds += 1
        vectors = r[self._j]-r[self._i]
        if self.period is not None:
            vectors[:, 2] += self._shifts*self.period
        radii = np.linalg.norm(vectors, axis=1)
        if np.any(radii <= 1e-10):
            raise ValueError('coincident atoms in SW configuration')
        active = radii < self.cutoff
        return self._i[active], self._j[active], vectors[active], radii[active]

    def evaluate(self, positions, *, representation='moments'):
        """Total energy and its Cartesian gradient, including fixed atoms."""
        r = np.asarray(positions, float)
        i, j, vectors, lengths = self._edges(r)
        if representation not in ('moments', 'direct'):
            raise ValueError('representation must be moments or direct')
        p, n = self.p, len(r)
        directions = vectors/lengths[:, None]
        gap = lengths-self.cutoff
        z = p['sigma']/gap
        expz = np.exp(z)
        power_p, power_q = (p['sigma']/lengths)**p['p'], (p['sigma']/lengths)**p['q']
        radial = p['B']*power_p-power_q
        pair = .5*p['A']*p['epsilon']*radial*expz
        pair_derivative = .5*p['A']*p['epsilon']*expz*(
            (-p['p']*p['B']*power_p+p['q']*power_q)/lengths
            - radial*p['sigma']/gap**2)
        edge_gradient = pair_derivative[:, None]*directions
        site_energy = np.bincount(i, weights=pair, minlength=n)
        w = np.exp(p['gamma']*z)
        c, strength = p['costheta0'], p['lambda']*p['epsilon']
        dw_dr = -p['gamma']*p['sigma']/gap**2*w
        if representation == 'moments':
            rho = np.bincount(i, weights=w, minlength=n)
            diagonal = np.bincount(i, weights=w*w, minlength=n)
            vector = np.stack([np.bincount(i, weights=w*directions[:, a], minlength=n)
                               for a in range(3)], axis=1)
            tensor = np.empty((n, 3, 3))
            for a in range(3):
                for b in range(3):
                    tensor[:, a, b] = np.bincount(i,
                        weights=w*directions[:, a]*directions[:, b], minlength=n)
            angular = (.5*np.einsum('iab,iab->i', tensor, tensor)
                       - c*np.einsum('ia,ia->i', vector, vector)
                       + .5*c*c*rho*rho-.5*(1-c)**2*diagonal)
            tn = np.einsum('iab,ib->ia', tensor[i], directions)
            vn = np.einsum('ia,ia->i', vector[i], directions)
            nn = np.einsum('ia,ia->i', tn, directions)
            d_energy_dw = nn-2*c*vn+c*c*rho[i]-(1-c)**2*w
            tangent = tn-c*vector[i]
            tangent -= directions*np.einsum('ia,ia->i', directions, tangent)[:, None]
            edge_gradient += strength*(
                (dw_dr*d_energy_dw)[:, None]*directions+2*(w/lengths)[:, None]*tangent)
            site_energy += strength*angular
        else:
            # Independent explicit centered unordered triples, no moment identity.
            counts = np.bincount(i, minlength=n)
            ends = np.cumsum(counts)
            starts = ends-counts
            for left in range(int(counts.max(initial=0))):
                for right in range(left+1, int(counts.max(initial=0))):
                    centers = np.flatnonzero(counts > right)
                    a, b = starts[centers]+left, starts[centers]+right
                    cosine = np.einsum('ij,ij->i', directions[a], directions[b])
                    delta = cosine-c
                    prefactor = strength*w[a]*w[b]
                    np.add.at(site_energy, centers, prefactor*delta**2)
                    for u, v in ((a, b), (b, a)):
                        derivative = prefactor[:, None]*(
                            (-p['gamma']*p['sigma']*delta**2/gap[u]**2)[:, None]*directions[u]
                            + (2*delta/lengths[u])[:, None]*(directions[v]-cosine[:, None]*directions[u]))
                        np.add.at(edge_gradient, u, derivative)
        gradient = np.stack([np.bincount(j, weights=edge_gradient[:, a], minlength=n)
            - np.bincount(i, weights=edge_gradient[:, a], minlength=n) for a in range(3)], axis=1)
        if not np.all(np.isfinite(gradient)) or not np.all(np.isfinite(site_energy)):
            raise FloatingPointError('nonfinite SW energy or gradient')
        return AtomicEvaluation(float(site_energy.sum()), gradient, site_energy)

    def hessian(self, positions, *, free_atoms=None):
        """Sparse analytic Hessian from independent per-site second-order jets.

        Fixed centers' energies are included before restricting the DOFs. This
        is a Cartesian Hessian, not a dynamical matrix or kinetic calibration.
        """
        from scipy.sparse import coo_matrix
        r = np.asarray(positions, float)
        i, j, vectors, _ = self._edges(r)
        free = np.arange(len(r)) if free_atoms is None else np.asarray(free_atoms, int)
        if (free.ndim != 1 or len(np.unique(free)) != len(free)
                or np.any(free < 0) or np.any(free >= len(r))):
            raise ValueError('distinct valid free atom indices required')
        mapping = np.full(len(r), -1, int)
        mapping[free] = np.arange(len(free))
        rows, cols, values = [], [], []
        starts = np.r_[0, np.cumsum(np.bincount(i, minlength=len(r)))]
        for center in range(len(r)):
            lo, hi = starts[center:center+2]
            neighbors = j[lo:hi]
            if not len(neighbors) or (mapping[center] < 0 and np.all(mapping[neighbors] < 0)):
                continue
            jet = environment_jet(vectors[lo:hi], self.p, 'direct')
            # Each neighbor-vector coordinate is R_neighbor - R_center.
            local, global_ids, signs = [], [], []
            for k, neighbor in enumerate(neighbors):
                for atom, sign in ((neighbor, 1.), (center, -1.)):
                    if mapping[atom] >= 0:
                        local.extend(range(3*k, 3*k+3))
                        global_ids.extend(range(3*mapping[atom], 3*mapping[atom]+3))
                        signs.extend([sign]*3)
            local, global_ids, signs = np.asarray(local), np.asarray(global_ids), np.asarray(signs)
            if not len(local):
                continue
            block = jet.hessian[np.ix_(local, local)]*np.outer(signs, signs)
            rows.extend(np.repeat(global_ids, len(global_ids)))
            cols.extend(np.tile(global_ids, len(global_ids)))
            values.extend(block.ravel())
        return coo_matrix((values, (rows, cols)), shape=(3*len(free),)*2).tocsr()


@dataclass
class CrackStrip:
    positions: np.ndarray
    fixed: np.ndarray
    front_groups: np.ndarray
    front_period: float
    basis_period: float
    width: float
    height: float
    crossing_bonds: np.ndarray
    frame: np.ndarray

    def seed(self, opening, *, tip=0., transition=6.):
        """Mode-I displacement seed; fixed grips are taken from this field.

        Far behind the tip two halves separate; ahead they stretch. No pair or
        triplet is removed. The smooth seed is not an elastic KI solution.
        """
        if (not np.all(np.isfinite([opening, tip, transition]))
                or opening < 0 or transition <= 0):
            raise ValueError('finite nonnegative opening and positive transition required')
        r = self.positions.copy()
        weight = .5*(1-np.tanh((r[:, 0]-tip)/transition))
        r[:, 1] += opening*((1-weight)*r[:, 1]/self.height+.5*weight*np.sign(r[:, 1]))
        return r

    def bond_gaps(self, positions):
        """Normal separations of the original cross-plane nearest neighbors.

        This observable can decrease on healing. A gap is not an irreversible
        fracture flag, probability, or rule for deleting an interaction.
        """
        bonds = self.crossing_bonds
        return np.asarray(positions)[bonds[:, 1], 1]-np.asarray(positions)[bonds[:, 0], 1]


def diamond_crack_strip(lattice, *, nx=8, ny=4, nz=4, grip_width=4.):
    """Diamond (111) shuffle strip, propagation [11-2], front [1-10]."""
    if (not np.isfinite(lattice) or lattice <= 0 or not np.isfinite(grip_width)
            or grip_width <= 0 or any(int(v) != v or v < 1 for v in (nx, ny, nz))):
        raise ValueError('positive lattice, grips and integer repeats required')
    nx, ny, nz = int(nx), int(ny), int(nz)
    frame = np.array([[1, 1, -2], [1, 1, 1], [1, -1, 0.]])
    frame /= np.linalg.norm(frame, axis=1)[:, None]
    cubic_cell = np.array([[.5, .5, -1], [1, 1, 1], [.5, -.5, 0.]])
    fcc = np.array([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
    diamond = np.concatenate([fcc, fcc+.25])
    translations = np.array(list(product(range(-3, 4), repeat=3)))
    fractional = (translations[:, None, :]+diamond).reshape(-1, 3)@np.linalg.inv(cubic_cell)
    basis = fractional[np.all((fractional >= -1e-10) & (fractional < 1-1e-10), axis=1)]
    if basis.shape != (12, 3):
        raise AssertionError('orthorhombic diamond cell must contain 12 atoms')
    basis[np.abs(basis) < 1e-10] = 0.
    copies = np.array(list(product(range(nx), range(ny), range(nz))))
    cell = lattice*cubic_cell@frame.T
    r = ((copies[:, None, :]+basis)@cell).reshape(-1, 3)
    width, height, period = np.diag(cell)*[nx, ny, nz]
    r[:, 0] -= .5*width
    heights = np.unique(np.round(r[:, 1], 9))
    spaces = np.diff(heights)
    shuffle = np.flatnonzero(spaces > .5*(spaces.min()+spaces.max()))
    cut = shuffle[np.argmin(abs((heights[shuffle]+heights[shuffle+1])*.5-height/2))]
    r[:, 1] -= .5*(heights[cut]+heights[cut+1])
    xmin, xmax, ymin, ymax = r[:, 0].min(), r[:, 0].max(), r[:, 1].min(), r[:, 1].max()
    fixed = ((r[:, 0] < xmin+grip_width) | (r[:, 0] > xmax-grip_width)
             | (r[:, 1] < ymin+grip_width) | (r[:, 1] > ymax-grip_width))
    if np.all(fixed):
        raise ValueError('strip has no mobile interior; increase repeats or reduce grips')
    groups = ((copies[:, 0]*ny+copies[:, 1])*len(basis))[:, None]+np.arange(len(basis))
    groups = np.broadcast_to(groups, (len(copies), len(basis))).ravel()
    # Shuffle bonds are normal to (111), with both ends in the same periodic image.
    pairs = cKDTree(r).query_pairs(lattice*np.sqrt(3)/4*1.01, output_type='ndarray')
    crossing = pairs[(r[pairs[:, 0], 1]*r[pairs[:, 1], 1] < 0)]
    swap = r[crossing[:, 0], 1] > 0
    crossing[swap] = crossing[swap, ::-1]
    return CrackStrip(r, fixed, groups, float(period), float(cell[2, 2]),
                      float(width), float(height), crossing, frame)


class RelaxedCoordinates:
    """Linear Cartesian reduction with optional coherent-front control.

    Independent displacements (default) allow local crack-front variation. The
    control ties equivalent atoms in different front repeats. Optional relative
    coordinates fix components of R_j-R_i exactly, leaving pair means free.
    All energies remain TOTAL energies of the actual periodic supercell.
    """
    def __init__(self, reference, fixed, *, groups=None, bond=None, components=(1,)):
        self.reference = np.asarray(reference, float).copy()
        fixed = np.asarray(fixed, bool)
        if self.reference.ndim != 2 or self.reference.shape[1] != 3 or fixed.shape != (len(reference),):
            raise ValueError('reference (n,3) and fixed (n,) required')
        raw = np.arange(len(reference)) if groups is None else np.asarray(groups)
        if raw.shape != fixed.shape or not np.all(np.isfinite(raw)):
            raise ValueError('one finite group label per atom required')
        _, self.labels = np.unique(raw, return_inverse=True)
        count = int(self.labels.max())+1
        group_fixed = np.zeros(count, bool)
        np.logical_or.at(group_fixed, self.labels, fixed)
        if np.any(group_fixed[self.labels] != fixed):
            raise ValueError('a coherent group cannot mix fixed and mobile atoms')
        self.mobile_groups = np.flatnonzero(~group_fixed)
        self.mapping = np.full(count, -1, int)
        self.mapping[self.mobile_groups] = np.arange(len(self.mobile_groups))
        self.nfree = len(self.mobile_groups)*3
        self.constraints = []
        if bond is not None:
            ids = np.asarray(bond, int)
            if ids.shape != (2,) or np.any(ids < 0) or np.any(ids >= len(reference)) or ids[0] == ids[1]:
                raise ValueError('two distinct valid bond atom indices required')
            g = self.mapping[self.labels[ids]]
            if np.any(g < 0) or g[0] == g[1]:
                raise ValueError('bond endpoints must be different mobile groups')
            if len(set(components)) != len(components) or any(c not in (0, 1, 2) for c in components):
                raise ValueError('distinct Cartesian constraint components required')
            for component in components:
                self.constraints.append((3*g[0]+component, 3*g[1]+component,
                    self.reference[ids[1], component]-self.reference[ids[0], component]))
        # The orthonormal difference coordinate is removed, its mean retained.
        self.keep = np.ones(self.nfree, bool)
        for _, b, _ in self.constraints:
            self.keep[b] = False
        self.dimension = int(self.keep.sum())

    def positions(self, variables, values=()):
        x, q = np.asarray(variables, float), np.asarray(values, float)
        if x.shape != (self.dimension,) or q.shape != (len(self.constraints),):
            raise ValueError('coordinate/constraint dimension mismatch')
        if not np.all(np.isfinite(x)) or not np.all(np.isfinite(q)):
            raise ValueError('finite coordinates required')
        displacement = np.zeros(self.nfree)
        displacement[self.keep] = x
        for target, (a, b, original) in zip(q, self.constraints):
            mean = displacement[a]/np.sqrt(2)
            displacement[a], displacement[b] = mean-.5*(target-original), mean+.5*(target-original)
        grouped = np.zeros((len(self.mapping), 3))
        grouped[self.mobile_groups] = displacement.reshape(-1, 3)
        return self.reference+grouped[self.labels]

    def pullback(self, gradient):
        grouped = np.zeros((len(self.mapping), 3))
        np.add.at(grouped, self.labels, gradient)
        flat = grouped[self.mobile_groups].ravel()
        reactions = []
        for a, b, _ in self.constraints:
            reactions.append(.5*(flat[b]-flat[a]))
            flat[a] = (flat[a]+flat[b])/np.sqrt(2)
        return flat[self.keep], np.asarray(reactions)

    def encode(self, positions):
        delta = np.asarray(positions)-self.reference
        grouped = np.zeros((len(self.mapping), 3))
        np.add.at(grouped, self.labels, delta)
        grouped /= np.bincount(self.labels)[:, None]
        flat = grouped[self.mobile_groups].ravel()
        for a, b, _ in self.constraints:
            flat[a] = (flat[a]+flat[b])/np.sqrt(2)
        return flat[self.keep]

    def tangent_matrix(self):
        """Constant dR/dx for the retained free coordinates (sparse)."""
        from scipy.sparse import coo_matrix
        lookup = np.full(self.nfree, -1, int)
        lookup[self.keep] = np.arange(self.dimension)
        mean_coordinates = {a: (a, b) for a, b, _ in self.constraints}
        mean_coordinates.update({b: (a, b) for a, b, _ in self.constraints})
        rows, cols, data = [], [], []
        for atom, label in enumerate(self.labels):
            group = self.mapping[label]
            if group < 0:
                continue
            for component in range(3):
                dof = 3*group+component
                if dof in mean_coordinates:
                    column, factor = lookup[mean_coordinates[dof][0]], 1/np.sqrt(2)
                else:
                    column, factor = lookup[dof], 1.
                rows.append(3*atom+component)
                cols.append(column)
                data.append(factor)
        return coo_matrix((data, (rows, cols)),
                          shape=(self.reference.size, self.dimension)).tocsr()


def relax_atoms(energy, coordinates, *, initial=None, values=(), tolerance=1e-7, maxiter=2500):
    """Force-checked static constrained relaxation, with raw optimizer status."""
    if not np.isfinite(tolerance) or tolerance <= 0 or int(maxiter) != maxiter or maxiter < 1:
        raise ValueError('positive finite force tolerance and iteration budget required')
    x = np.zeros(coordinates.dimension) if initial is None else coordinates.encode(initial)
    sites_zero = energy.evaluate(coordinates.positions(x, values)).site_energy.copy()

    def objective(variables):
        current = energy.evaluate(coordinates.positions(variables, values))
        gradient, _ = coordinates.pullback(current.gradient)
        # Subtract per site before summing: avoid subtracting two extensive
        # totals when a tiny relaxation step is being resolved.
        return float(np.sum(current.site_energy-sites_zero)), gradient

    result = minimize(objective, x, jac=True, method='L-BFGS-B',
        options={'gtol':tolerance*.25, 'ftol':0., 'maxiter':int(maxiter), 'maxls':40, 'maxcor':20})
    positions = coordinates.positions(result.x, values)
    evaluation = energy.evaluate(positions)
    projected, reactions = coordinates.pullback(evaluation.gradient)
    residual = float(np.max(np.abs(projected), initial=0))
    return positions, {
        'energy_eV':evaluation.energy, 'free_gradient_max_eV_A':residual,
        'force_tolerance_eV_A':tolerance, 'force_converged':bool(residual <= tolerance),
        'constraint_reaction_eV_A':reactions.tolist(),
        'optimizer_success':bool(result.success), 'optimizer_message':str(result.message),
        'iterations':int(result.nit), 'evaluations':int(result.nfev),
        'stability_certified':False, 'physical_time_seconds':None,
    }
