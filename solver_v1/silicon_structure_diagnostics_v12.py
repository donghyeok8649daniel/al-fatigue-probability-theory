"""Geometric evidence for intact Si configurations, never a crack classifier.

Reference-neighbor local affine residuals, graph bottlenecks and sectional pair
counts use explicit distance conventions. Results do not define physical bonds,
plastic strain, an initiation boundary, or statistical transition rates.
"""
from __future__ import annotations
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import maximum_flow, breadth_first_order


def pair_set(positions,cutoff):
    r=np.asarray(positions,float)
    if r.ndim!=2 or r.shape[1]!=3 or not np.all(np.isfinite(r)) or not np.isfinite(cutoff) or cutoff<=0:
        raise ValueError('finite 3D coordinates and positive cutoff required')
    return {tuple(x) for x in cKDTree(r).query_pairs(float(cutoff),output_type='ndarray').tolist()}


def neighbor_lists(n,pairs):
    neighbors=[[] for _ in range(n)]
    for i,j in sorted(pairs):
        if not 0<=i<j<n:raise ValueError('unique ascending atom pairs required')
        neighbors[i].append(j);neighbors[j].append(i)
    return neighbors


def local_affine(before,after,neighbors):
    """Least squares relative-neighbor map X A = Y. No rank-deficient strain.

    D2 is the mean squared residual per neighbor, with Angstrom^2 units.
    Green strain eigenvalues come from (A A^T - I)/2 for row-vector coordinates.
    Initial neighborhoods and central atom labels are kept fixed.
    """
    x=np.asarray(before,float);y=np.asarray(after,float)
    if x.shape!=y.shape or x.ndim!=2 or x.shape[1]!=3 or len(neighbors)!=len(x):
        raise ValueError('matching coordinates and neighborhoods required')
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):raise ValueError('nonfinite coordinate')
    records=[]
    for i,indices in enumerate(neighbors):
        ids=np.asarray(indices,int)
        if len(ids)==0:
            records.append(dict(atom=i,neighbors=0,rank=0,residual_degrees_of_freedom=0,condition=None,D2_A2=None,green_principal_strains=None));continue
        dx=x[ids]-x[i];dy=y[ids]-y[i]
        matrix,_,rank,singular=np.linalg.lstsq(dx,dy,rcond=None)
        residual=dy-dx@matrix
        valid=rank==3
        records.append(dict(atom=i,neighbors=len(ids),rank=int(rank),residual_degrees_of_freedom=int(3*(len(ids)-rank)),
            condition=float(singular[0]/singular[-1]) if valid else None,
            D2_A2=float(np.mean(np.sum(residual**2,axis=1))),
            green_principal_strains=np.linalg.eigvalsh((matrix@matrix.T-np.eye(3))/2).tolist() if valid else None))
    return records


def edge_bottleneck(n,pairs,lower,upper):
    """Number of edge-disjoint cutoff-graph paths between two rigid grips.

    Unit capacities are topological only. A min cut is NOT a fracture path,
    traction law, measured strength or bond-energy sum.
    """
    lower=np.asarray(lower,bool);upper=np.asarray(upper,bool)
    if lower.shape!=(n,) or upper.shape!=(n,) or np.any(lower&upper) or not lower.any() or not upper.any():
        raise ValueError('disjoint nonempty masks required')
    row=[];col=[];data=[];large=2*len(pairs)+1
    for i,j in sorted(pairs):
        row.extend([i,j]);col.extend([j,i]);data.extend([1,1])
    for i in np.flatnonzero(lower):row.append(n);col.append(i);data.append(large)
    for i in np.flatnonzero(upper):row.append(i);col.append(n+1);data.append(large)
    capacity=coo_matrix((np.asarray(data,dtype=np.int64),(row,col)),shape=(n+2,n+2)).tocsr()
    flow=maximum_flow(capacity,n,n+1)
    residual=(capacity-flow.flow).tocsr();residual.data[residual.data<=0]=0;residual.eliminate_zeros()
    reached=breadth_first_order(residual,n,directed=True,return_predecessors=False)
    source_side=np.zeros(n+2,bool);source_side[reached]=True
    cut=[(i,j) for i,j in sorted(pairs) if source_side[i]!=source_side[j]]
    if len(cut)!=int(flow.flow_value):raise ValueError('max-flow/min-cut identity failed')
    return dict(edge_disjoint_paths=int(flow.flow_value),cut_pairs=cut,
        source_side_atoms=np.flatnonzero(source_side[:n]).tolist(),cut_is_physical_crack=False)


def reference_sections(reference,current,initial_pairs,current_pairs,free):
    """Cuts between reference z layers, keeping reference atom membership.

    Counts are reference-labelled sections. Their physical geometry can warp;
    counts and mean plane separation are not a local stress/area estimate.
    """
    r=np.asarray(reference,float);q=np.asarray(current,float);free=np.asarray(free,bool)
    # Rounding groups exactly intended crystal planes despite floating point noise.
    levels=np.unique(np.round(r[:,2],10));rows=[]
    for za,zb in zip(levels[:-1],levels[1:]):
        z=(za+zb)/2;side=r[:,2]<z
        if not np.any(free&side) or not np.any(free&~side):continue
        old={p for p in initial_pairs if side[p[0]]!=side[p[1]]}
        new={p for p in current_pairs if side[p[0]]!=side[p[1]]}
        layer_a=np.isclose(r[:,2],za,rtol=0,atol=1e-8)
        layer_b=np.isclose(r[:,2],zb,rtol=0,atol=1e-8)
        rows.append(dict(reference_midpoint_A=float(z),reference_gap_A=float(zb-za),
            current_mean_layer_gap_A=float(q[layer_b,2].mean()-q[layer_a,2].mean()),
            initial_crossing_pairs=len(old),current_crossing_pairs=len(new),
            lost_initial_pairs=len(old-new),new_pairs=len(new-old)))
    return rows
