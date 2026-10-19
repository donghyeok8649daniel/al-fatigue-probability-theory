"""Run isolated analytical/numerical demonstrations for Parts I--III.

No solver settings, material calibration or historical results are changed.
Every example has dimensionless illustrative parameters. No fatigue-life
Monte Carlo, aluminium prediction or physical time calibration is performed.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import sys
import json
import time
import traceback
import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
from scipy.linalg import block_diag, solve_continuous_lyapunov
from foundations import (B, SIGMA, row_direct, row_reciprocal, row_tail_value_bound,
    lj_row, ZigzagCut, plane_direct, plane_reciprocal, plane_tail_bound,
    finite_crossing_counts, finite_crossing_formula)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=False)
    records=[];started=time.monotonic()
    def check(name,fn):
        t=time.monotonic()
        try:
            details=fn();passed=True;error=None
        except Exception as exc:
            details={};passed=False;error=repr(exc)
            (out/(name+'_failure.txt')).write_text(traceback.format_exc(),encoding='utf-8')
        records.append(dict(name=name,passed=passed,seconds=time.monotonic()-t,error=error,details=details))
        print(('PASS ' if passed else 'FAIL ')+name,flush=True)

    def rows():
        result=[]
        for a,s in ((.12,.49),(.45,.13),(.8,-.21),(1.6,.31),(5.,-.07)):
            for p in (3,6):
                exact=row_reciprocal(a,s,p,modes=192)
                reference=row_direct(a,s,p,radius=16384)
                scale=np.maximum(1.,np.abs(reference))
                error=np.abs(exact-reference)/scale
                assert np.max(error)<2.e-12,(a,s,p,error)
                low=row_reciprocal(a,s,p,modes=96)
                assert np.max(np.abs(exact-low)/scale)<2.e-12
                result.append(dict(a=a,s=s,p=p,maximumScaledError=float(np.max(error)),
                    valueAbsoluteError=float(abs(exact[0]-reference[0])),realValueTailBound=row_tail_value_bound(p,16384),
                    realRadius=16384,reciprocalModes=192))
        # This direct sum is independent of production Fourier code; the one
        # source comparison only locates the established model's meaning.
        root=Path(__file__).resolve().parents[2];sys.path.insert(0,str(root))
        from solver_v1.lattice_bessel import two_row_lj_infinite_energy_gradient,FourierLatticeConfig
        own=lj_row(.81,.16,modes=96)
        prod=two_row_lj_infinite_energy_gradient(.81,.16,epsilon=1.,sigma_lj=SIGMA,b=B,
                                               config=FourierLatticeConfig(tol=1.e-14,max_modes=96))
        assert np.allclose(own[:3],prod[:3],rtol=2.e-12,atol=2.e-12)
        return dict(points=result,productionComparison=dict(model='original one-upper-site/infinite-lower-row only',
            energyGradientMaximumDifference=float(np.max(np.abs(own[:3]-np.asarray(prod[:3]))))),
            forceSign='Internal force is minus reported gradient; Hessian channels checked separately from force.')
    check('row_energy_force_hessian_direct_vs_poisson',rows)

    def planes():
        records=[]
        for p in (3,6):
            for d,delta in ((.55,(.13,-.11)),(1.3,(.25,.2))):
                reciprocal=plane_reciprocal(d,delta,p,shell=16)
                convergence=[]
                for radius in (64,128,256):
                    real=plane_direct(d,delta,p,shell=radius)
                    value_err=abs(real[0]-reciprocal[0]);bound=plane_tail_bound(p,radius,delta)
                    all_error=max(value_err,float(np.max(abs(real[1]-reciprocal[1]))),float(np.max(abs(real[2]-reciprocal[2]))))
                    assert value_err<=bound+3.e-12*max(1.,abs(reciprocal[0]))
                    assert all_error<max(2.e-10*max(1.,abs(reciprocal[0])),30*bound)
                    convergence.append(dict(realShell=radius,valueError=value_err,valueTailBound=bound,
                                            largestValueGradientHessianError=all_error))
                higher=plane_reciprocal(d,delta,p,shell=24)
                assert abs(higher[0]-reciprocal[0])<1.e-11*max(1.,abs(higher[0]))
                records.append(dict(p=p,d=d,delta=delta,order='K_'+str(p-1),reciprocalShell=16,convergence=convergence))
        return dict(records=records,normalization='Per triangular primitive plane cell; area=sqrt(3)b^2/2, reciprocal prefactor1/area.')
    check('planar_kernel_order_normalization_and_real_tail_refinement',planes)

    # Find a stress-free fixed-b reference for this declared row stack only.
    h=brentq(lambda gap:ZigzagCut(gap).features(0.,0.)[1],.45,1.3,xtol=2.e-14)
    cut=ZigzagCut(h,layers=16,modes=64)
    pristine=cut.features(0.,0.);H=np.array([[pristine[3],pristine[4]],[pristine[4],pristine[5]]])
    assert np.min(np.linalg.eigvalsh(H))>0

    def counting():
        L,U=4,5;counts=finite_crossing_counts(L,U)
        assert sum(counts.values())==L*U
        assert all(v==finite_crossing_formula(k,L,U) for k,v in counts.items())
        assert [counts[k] for k in range(1,9)]==[1,2,3,4,4,3,2,1]
        for k in range(1,13):
            assert finite_crossing_counts(16,16)[k]==k
        delta,s=.09,.13;direct=0.;fourier=0.
        for low in range(1-L,1):
            for high in range(1,U+1):
                k=high-low;seff=s+(1-k)*B/2;ref_s=(1-k)*B/2
                def raw(d,slip):
                    return 4*(SIGMA**12*row_direct(d,slip,6,radius=2048)[0]-SIGMA**6*row_direct(d,slip,3,radius=2048)[0])
                direct+=raw(k*h+delta,seff)-raw(k*h,ref_s)
        for k,count in counts.items():
            fourier+=count*(lj_row(k*h+delta,s+(1-k)*B/2)[0]-lj_row(k*h,(1-k)*B/2)[0])
        assert abs(direct-fourier)<2.e-11
        # An incorrect extra chain-rule k is measurably different from the
        # derivative of an active cut. Homogeneous expansion is another path.
        true_force=sum(count*lj_row(k*h+delta,s+(1-k)*B/2)[1] for k,count in counts.items())
        wrong_force=sum(count*k*lj_row(k*h+delta,s+(1-k)*B/2)[1] for k,count in counts.items())
        assert abs(true_force-wrong_force)>1.e-3
        return dict(lowerRows=L,upperRows=U,counts=counts,crossingPairs=L*U,directEnergy=direct,
                    countedBesselEnergy=fourier,absoluteDifference=abs(direct-fourier),
                    cutDerivative=true_force,incorrectHomogeneousChainRule=wrong_force)
    check('explicit_cut_bond_enumeration_and_path_derivative',counting)

    def cleavage():
        assert abs(pristine[0])<1.e-13 and np.max(abs(pristine[1:3]))<2.e-10
        q=(.08,.14);base=cut.features(*q);refined=ZigzagCut(h,layers=24,modes=96).features(*q)
        assert np.max(abs(base-refined))<2.e-10
        eps=2.e-5
        grad=np.array([(cut.features(q[0]+eps,q[1])[0]-cut.features(q[0]-eps,q[1])[0])/(2*eps),
                       (cut.features(q[0],q[1]+eps)[0]-cut.features(q[0],q[1]-eps)[0])/(2*eps)])
        hess_a=(cut.features(q[0]+eps,q[1])[1:3]-cut.features(q[0]-eps,q[1])[1:3])/(2*eps)
        hess_s=(cut.features(q[0],q[1]+eps)[1:3]-cut.features(q[0],q[1]-eps)[1:3])/(2*eps)
        expected=np.array([[base[3],base[4]],[base[4],base[5]]])
        assert np.max(abs(grad-base[1:3]))<1.e-5
        assert np.max(abs(np.column_stack((hess_a,hess_s))-expected))<2.e-4
        endpoint=cut.separation_limit();approach=[]
        for opening in (4.,8.,16.,32.,64.):
            energy=cut.features(opening,0.)[0]
            approach.append(dict(opening=opening,energy=energy,endpointError=abs(endpoint-energy)))
        assert all(a['endpointError']>b['endpointError'] for a,b in zip(approach,approach[1:]))
        return dict(geometry='2D zigzag half-stacks; per lateral repeat, not J/m²',h=h,
            referenceFeatures=pristine.tolist(),positiveHessianEigenvalues=np.linalg.eigvalsh(H).tolist(),
            layerModeRefinementMaximumDifference=float(np.max(abs(base-refined))),
            finiteDifferenceGradientError=float(np.max(abs(grad-base[1:3]))),
            finiteDifferenceHessianError=float(np.max(abs(np.column_stack((hess_a,hess_s))-expected))),
            separationLimit=endpoint,approach=approach,zeroMode='Infinite weighted Hurwitz-zeta sum retained exactly.')
    check('counted_rigid_cut_derivatives_refinement_and_separation_limit',cleavage)

    def pmf():
        T=.2
        def stiffness(q):return 2*np.exp(.8*q[0]+.4*(1-np.cos(2*np.pi*q[1])))
        q0=(0.,0.);z0=quad(lambda y:np.exp(-stiffness(q0)*y*y/(2*T)),-np.inf,np.inf,epsabs=1.e-12)[0]
        rows=[]
        for q in ((.1,.17),(-.03,.08),(.2,-.21)):
            k=stiffness(q);z,err=quad(lambda y:np.exp(-k*y*y/(2*T)),-np.inf,np.inf,epsabs=1.e-12)
            delta_w=cut.features(*q)[0];entropy=-T*np.log(z/z0);analytic=T/2*np.log(k/stiffness(q0))
            assert abs(entropy-analytic)<2.e-12 and abs(z-np.sqrt(2*np.pi*T/k))<2.e-12
            rows.append(dict(q=q,rigidAndRelaxedEnergy=delta_w,PMF=delta_w+entropy,
                             entropyCorrection=entropy,independentGaussianIntegralError=err))
        assert max(abs(row['entropyCorrection']) for row in rows)>.01
        return dict(T=T,model='V(q,y)=W_Gamma_Bessel(q)+kappa(q)y²/2',sameRigidAndRelaxedEnergy=True,
                    PMFEqualsRigid=False,records=rows,meaning='Exact analytical counterexample to universal W=A, not an Al entropy estimate.')
    check('same_bessel_energy_different_PMF_exact_entropy_counterexample',pmf)

    def finite_bath():
        masses=np.array([1.,1.5,.8]);frequencies=np.array([3.,5.,9.]);coupling=np.array([.7,-.4,.6]);mu=.4
        initial_q=.006;initial_p=.01;r0=np.array([.03,-.02,.01]);P0=np.array([.02,.03,-.01])
        x0=r0+coupling*initial_q/(masses*frequencies**2)
        initial=np.r_[initial_q,initial_p,x0,P0]
        def forcing(t):return .02*np.sin(2*t) # Normal external force only.
        def rhs(t,state):
            q,p=state[:2];x=state[2:5];P=state[5:]
            residual=x-coupling*q/(masses*frequencies**2)
            return np.r_[p/mu,-cut.features(q,0.)[1]+coupling@residual+forcing(t),P/masses,-masses*frequencies**2*residual]
        def kernel(t):return np.sum(coupling**2/(masses*frequencies**2)*np.cos(frequencies*t))
        def noise(t):return coupling@(r0*np.cos(frequencies*t)+P0/(masses*frequencies)*np.sin(frequencies*t))
        def internal_energy(state):
            q,p=state[:2];x=state[2:5];P=state[5:]
            return p*p/(2*mu)+cut.features(q,0.)[0]+np.sum(P*P/(2*masses)+masses*frequencies**2/2*(x-coupling*q/(masses*frequencies**2))**2)
        convergence=[]
        for tolerance,step in ((2.e-8,.02),(2.e-10,.01)):
            sol=solve_ivp(rhs,(0.,1.2),initial,method='DOP853',rtol=tolerance,atol=tolerance*.01,max_step=step,dense_output=True)
            assert sol.success
            residuals=[]
            for t in np.linspace(.06,1.2,14):
                state=sol.sol(t);conv=quad(lambda u:kernel(t-u)*sol.sol(u)[1]/mu,0,t,epsabs=2.e-11,limit=100)[0]
                gle=-cut.features(state[0],0.)[1]-conv+noise(t)+forcing(t)
                residuals.append(abs(gle-rhs(t,state)[1]))
            work=quad(lambda t:forcing(t)*sol.sol(t)[1]/mu,0,1.2,epsabs=2.e-11,limit=100)[0]
            work_error=abs(internal_energy(sol.y[:,-1])-internal_energy(initial)-work)
            assert max(residuals)<1.e-7 and work_error<1.e-8
            convergence.append(dict(rtol=tolerance,maximumStep=step,RHSEvaluations=sol.nfev,
                maximumHamiltonianVsGLEForceError=max(residuals),normalCyclicWorkBalanceError=work_error))
        T=.2;covariance_error=[]
        for t,u in ((.1,.7),(.3,.31),(1.,.2)):
            covariance=np.sum(coupling**2*(T/(masses*frequencies**2)*np.cos(frequencies*t)*np.cos(frequencies*u)
                +masses*T/(masses*frequencies)**2*np.sin(frequencies*t)*np.sin(frequencies*u)))
            covariance_error.append(abs(covariance-T*kernel(t-u)))
        recurrence=abs(kernel(2*np.pi)-kernel(0.));assert recurrence<1.e-14
        slip_error=abs(np.sum(-coupling**2*initial_q/(masses*frequencies**2)*np.cos(frequencies*.4))+kernel(.4)*initial_q)
        assert max(covariance_error)<1.e-14 and slip_error<1.e-14
        return dict(potential='Full nonlinear counted Bessel cut, s constrained to0; not its harmonic replacement',
            dynamics='Deterministic finite Hamiltonian bath under normal cyclic force; no fatigue trajectories',
            refinements=convergence,analyticCanonicalNoiseCovarianceError=max(covariance_error),
            GaussianButNotWhite='Finite Gaussian bath has exact 2pi kernel recurrence',recurrenceError=recurrence,
            unshiftedPreparationInitialSlipError=slip_error,materialParameters='illustrative dimensionless only')
    check('finite_hamiltonian_bath_exact_memory_FDT_and_initial_preparation',finite_bath)

    def controlled_limits():
        T=.07;mu=.3;eps=.1;friction=np.array([[2.,.3],[.3,1.4]]);I=np.eye(2);Z=np.zeros((2,2))
        drift=np.block([[Z,I/mu,Z],[-H,Z,-I],[Z,friction/(eps*mu),-I/eps]])
        noise=np.zeros((6,2));noise[4:]=np.sqrt(2*T)*np.linalg.cholesky(friction)/eps
        cov=solve_continuous_lyapunov(drift,-noise@noise.T)
        exact=block_diag(T*np.linalg.inv(H),T*mu*I,T*friction/eps)
        assert np.max(abs(cov-exact))<2.e-10
        assert np.max(np.real(np.linalg.eigvals(drift)))<0
        omega=1.;drive=np.array([1.,0.])
        def response(m,e):return np.linalg.solve(H-m*omega**2*I+1j*omega*friction/(1+1j*omega*e),drive)
        white=response(mu,0.);over=response(0.,0.)
        memory_errors=[float(np.linalg.norm(response(mu,e)-white)/np.linalg.norm(white)) for e in (.25,.05,.01)]
        inertial_errors=[float(np.linalg.norm(response(m,0.)-over)/np.linalg.norm(over)) for m in (.3,.06,.012)]
        assert all(a>b for a,b in zip(memory_errors,memory_errors[1:]))
        assert all(a>b for a,b in zip(inertial_errors,inertial_errors[1:]))
        mobility=np.linalg.inv(friction);rates=np.linalg.eigvals(mobility@H)
        assert np.allclose(np.linalg.eigvals((2*mobility)@H),2*rates)
        return dict(reference='Local quadratic expansion using the checked Bessel-cut Hessian; approximation explicitly local',
            Hessian=H.tolist(),friction=friction.tolist(),coloredEmbeddingEquilibriumCovarianceError=float(np.max(abs(cov-exact))),
            memoryCorrelationTimes=[.25,.05,.01],relativeColoredToWhiteResponseErrors=memory_errors,
            masses=[.3,.06,.012],relativeInertialToOverdampedResponseErrors=inertial_errors,
            samePMFDifferentRelaxationRates=dict(rates=rates.real.tolist(),doubledMobilityRates=(2*rates.real).tolist()),
            physicalInference='No measured scale separation, mobility or seconds inferred from this mathematical control.')
    check('memory_markov_and_inertia_limits_with_unchanged_equilibrium_PMF',controlled_limits)

    def multiplicative():
        T=.25;s=np.linspace(-.5,.5,513)
        features=np.array([cut.features(0.,v) for v in s]);U=features[:,0];Us=features[:,2]
        P=np.exp(-(U-np.min(U))/T);P/=np.trapezoid(P,s);Ps=-Us*P/T
        mobility=1.2+.25*np.sin(2*np.pi*s);dm=.5*np.pi*np.cos(2*np.pi*s)
        correct=(-mobility*Us+T*dm)*P-T*(dm*P+mobility*Ps)
        wrong=(-mobility*Us)*P-T*(dm*P+mobility*Ps)
        assert np.max(abs(correct))<2.e-12 and np.max(abs(wrong))>.01
        return dict(reference='Periodic registry quotient at fixed a; not a normalizable Gibbs law on unwrapped s',
            maximumCorrectGibbsFlux=float(np.max(abs(correct))),maximumMissingItoDriftFlux=float(np.max(abs(wrong))),
            rule='b_Ito=-M A_prime+kBT M_prime; all force comes from the same Bessel cut')
    check('multiplicative_mobility_drift_cannot_be_omitted',multiplicative)

    result=dict(scope='PART I-III independent mathematical demonstrations; not calibrated Al or production adoption',
        versions=dict(python=sys.version.split()[0],numpy=np.__version__,scipy=scipy.__version__),
        baseCommit='c8ccff7553cbb2e66b0c37a76ee8a24808478824',parameters=dict(epsilon=1.,sigma_lj=SIGMA,b=B,h=h),
        checks=records,passed=sum(r['passed'] for r in records),failed=sum(not r['passed'] for r in records),
        totalSeconds=time.monotonic()-started,productionSolverChanged=False,physicalTimeIdentified=False,
        gaussianBathDerivedForActualAl=False,fullExistingSolverOrAppSuiteRun=False)
    (out/'verification.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'],failed=result['failed'],seconds=result['totalSeconds'])),flush=True)
    return 1 if result['failed'] else 0


if __name__=='__main__':raise SystemExit(main())
