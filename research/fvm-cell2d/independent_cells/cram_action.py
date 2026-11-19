"""IPF CRAM48 sparse matrix exponential for stiff 2D FVM generators.

Coefficients: OpenMC deplete/cram.py, accessed 2026-11-18.
https://docs.openmc.org/en/stable/_modules/openmc/deplete/cram.html
Method: M. Pusa, Nucl. Sci. Eng. 182 (2016) 297-318, doi:10.13182/NSE15-26.
OpenMC coefficients are distributed under the MIT license in OPENMC_LICENSE.txt.
No OpenMC package dependency is needed.
"""
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu

POLES=np.array([
 -4.465731934165702e1,-5.284616241568964,-8.867715667624458,3.493013124279215,
 15.64102508858634,17.42097597385893,-28.34466755180654,16.61569367939544,
 8.011836167974721,-2.056267541998229,14.49208170441839,18.53807176907916,
 9.932562704505182,-22.44223871767187,.8590014121680897,-12.86192925744479,
 11.64596909542055,18.06076684783089,5.870672154659249,-35.42938819659747,
 19.01323489060250,18.85508331552577,-17.34689708174982,13.16284237125190
])+1j*np.array([
 62.33225190695437,40.57499381311059,43.25515754166724,32.81615453173585,
 15.58061616372237,10.76629305714420,54.92841024648724,13.16994930024688,
 27.80232111309410,37.94824788914354,17.99988210051809,5.974332563100539,
 25.32823409972962,51.79633600312162,35.36456194294350,46.00304902833652,
 22.87153304140217,8.368200580099821,30.29700159040121,58.34381701800013,
 1.194282058271408,3.583428564427879,48.83941101108207,20.42951874827759])
RESIDUES=np.array([
 638.7380733878774,190.9896179065730,423.6195226571914,464.5770595258726,
 776.5163276752433,1907.115136768522,2909.892685603256,194.4772206620450,
 138279.9786972332,5628.442079602433,215.1681283794220,1324.720240514420,
 16175.48476343347,111.2729040439685,107.4624783191125,88.35727765158191,
 93.54078136054179,94.18142823531573,104.0012390717851,68.61882624343235,
 87.66654491283722,105.6007619389650,77.38987569039419,104.1366366475571
])+1j*np.array([
 -674.3912502859256,-397.3203432721332,-2041.233768918671,-1652.917287299683,
 -17836.17639907328,-58870.68595142284,-9953.255345514560,-1427.131226068449,
 -3256885.197214938,-29242.84515884309,-1121.774011188224,-63700.88443140973,
 -1008798.413156542,-88.37109731680418,-145.7246116408180,-63.88286188419360,
 -219.5424319460237,-671.9055740098035,-169.3747595553868,-11.77598523430493,
 -4596.464999363902,-1738.294585524067,-43.11715386228984,-277.7743732451969])
ALPHA0=2.258038182743983e-47


def action(matrix,initial,duration,substeps=1):
    scaled=matrix.tocsc()*(duration/substeps)
    identity=sparse.eye(matrix.shape[0],format='csc')
    y=np.array(initial,dtype=float,copy=True)
    max_lu_bytes=0;max_solve_residual=0.
    # Retain one LU at a time. Apply it to every substep through stored stage
    # vectors only when substeps=1; substep runs below use separate full actions.
    if substeps>1:
        checks=[]
        for _ in range(substeps):
            y,c=action(matrix,y,duration/substeps,1);checks.append(c)
        return y,dict(maximum_LU_bytes=max(c['maximum_LU_bytes'] for c in checks),
            maximum_scaled_solve_residual=max(c['maximum_scaled_solve_residual'] for c in checks),
            rational_order=48,substeps=substeps)
    for alpha,pole in zip(RESIDUES,POLES):
        shifted=scaled-pole*identity
        lu=splu(shifted)
        rhs=y.astype(complex)
        z=lu.solve(rhs)
        error=abs(shifted@z-rhs).max()/max(1.,abs(rhs).max())
        max_solve_residual=max(max_solve_residual,float(error))
        max_lu_bytes=max(max_lu_bytes,sum(a.nbytes for mat in (lu.L,lu.U)
            for a in (mat.data,mat.indices,mat.indptr)))
        y+=2*np.real(alpha*z)
        del lu
    y*=ALPHA0
    return y,dict(maximum_LU_bytes=max_lu_bytes,maximum_scaled_solve_residual=max_solve_residual,
        rational_order=48,substeps=1)
