"""Two successive rational orders on the identical local SG generator.

IPF CRAM16/48 coefficients and MIT attribution: see cram_action.py and
OPENMC_LICENSE.txt. CRAM16 is only an error comparison, never the saved result.
"""
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu
import fvm2d as m
from cram_action import action

ALPHA16=np.array([
 5.464930576870210e3-3.797983575308356e4j,
 9.045112476907548e1-1.115537522430261e3j,
 2.344818070467641e2-4.228020157070496e2j,
 9.453304067358312e1-2.951294291446048e2j,
 7.283792954673409e2-1.205646080220011e5j,
 3.648229059594851e1-1.155509621409682e2j,
 2.547321630156819e1-2.639500283021502e1j,
 2.394538338734709e1-5.650522971778156j])
POLE16=np.array([
 3.509103608414918+8.436198985884374j,
 5.948152268951177+3.587457362018322j,
 -5.264971343442647+16.22022147316793j,
 1.419375897185666+10.92536348449672j,
 6.416177699099435+1.194122393370139j,
 4.993174737717997+5.996881713603942j,
 -1.413928462488886+13.49772569889275j,
 -10.84391707869699+19.27744616718165j])


def coarse_action(Q,initial,duration):
    identity=sparse.eye(Q.shape[0],format='csc');scaled=Q.tocsc()*duration
    y=np.array(initial,copy=True,dtype=float)
    for alpha,pole in zip(ALPHA16,POLE16):
        solver=splu(scaled-pole*identity)
        y+=2*np.real(alpha*solver.solve(y.astype(complex)))
        del solver
    return y*2.124853710495224e-16


def advance(operator,mass,cdf,duration):
    # The positive full-domain action is cheaper for moderately stiff phases.
    nu=1.01*float(np.max(-operator.matrix.diagonal()))
    if nu*duration<15_000:
        return m.advance_uniform_batch(operator,mass,cdf,duration)
    n=len(operator.keep)
    Q=sparse.bmat([[operator.matrix,sparse.csr_matrix((n,1))],
        [sparse.csr_matrix(operator.event[None,:]),sparse.csr_matrix((1,1))]],format='csc')
    mask=np.ones(operator.common_states,dtype=bool);mask[operator.keep]=False
    moving=mass[:,mask].sum(axis=1)
    initial=np.vstack([mass[:,operator.keep].T,np.zeros((1,len(mass)))])
    coarse=coarse_action(Q,initial,duration)
    fine,check=action(Q,initial,duration,1)
    result=np.zeros_like(mass);result[:,operator.keep]=fine[:-1].T
    final_f=cdf+moving+fine[-1]
    l1=float(abs(fine-coarse).sum(axis=0).max());df=float(abs(fine[-1]-coarse[-1]).max())
    balance=float(abs(result.sum(axis=1)+final_f-1.).max())
    negative=float((-np.minimum(result,0.)).sum(axis=1).max())
    if not (l1<2e-12 and df<2e-15 and balance<1e-9 and negative<2e-12):
        raise ArithmeticError(dict(method='CRAM16_vs_48',L1=l1,CDF=df,balance=balance,negative=negative))
    if not np.all(final_f>=cdf-2e-15):raise ArithmeticError('Nonmonotone first passage')
    return result,final_f,dict(method='full_sparse_CRAM48_successive_orders_16_48',state_dimension=2,
        independent_cell_count=len(mass),mass_table_shape=list(mass.shape),
        common_operator_shape=list(operator.matrix.shape),maximum_successive_L1=l1,
        maximum_successive_CDF=df,maximum_mass_residual=balance,maximum_negative_mass=negative,
        explicit_moving_event_increment=moving.tolist(),first_passage_flux=fine[-1].tolist(),
        numerical_mass_loss_added_to_event=False,probability_renormalized=False,
        probability_table_bytes=result.nbytes,rational_action_checks=[check])
