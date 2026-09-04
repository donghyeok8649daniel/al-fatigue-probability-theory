# Candidate physical construction of the initial spacing probability P0

Status: **CANDIDATE SUPPORTING THEORY — for the local-traction P0 -> P propagator; not yet promoted to the active exact finite-chain theory.**

Purpose: construct the initial spacing density from a physical specimen state rather than selecting a named probability family.

The hard reduced-model requirement is

$$
P_0(a\mid x)+\sigma_n(x,0:t)\longrightarrow P_a(a,t\mid x).
$$

A central consistency issue is that an instantaneous thermal spacing distribution and a static structural spacing distribution are not the same object. They must not be mixed.

## 1. Recommended interpretation for the strict P0-only candidate

For the strict P0-only local-traction candidate, define P0 as a **slow structural / coarse-grained inter-layer-spacing distribution** after averaging over fast phonon motion during specimen preparation.

Denote this candidate structural density by

$$
P_0^{\mathrm{str}}(a).
$$

Its compatible preparation condition is that the coarse-grained local spacing rates vanish at the chosen initial time,

$$
c_0=0.
$$

Then

$$
F_0^{\mathrm{str}}(a,c)=P_0^{\mathrm{str}}(a)\,\delta(c),
$$

which preserves the desired P0-only initialization of the candidate propagator.

This definition does not claim that real atoms have zero thermal velocity. It defines the reduced state after the fast thermal motion has been averaged out.

## 2. Exact spatial push-forward from a measured or computed structural spacing field

Let the specimen region represented by one probability state be Omega and let w(x) be a normalized sampling weight,

$$
\int_{\Omega}w(x)\,dx=1.
$$

Let the prepared local structural spacing be a_0^{\mathrm{str}}(x). Then the initial density is the push-forward

$$
P_0^{\mathrm{str}}(a)
=\int_{\Omega}w(x)\,\delta[a-a_0^{\mathrm{str}}(x)]\,dx.
$$

No Gaussian, Weibull, Boltzmann, or other named density family is introduced. If the measured structural spacing field is known at discrete sample points x_j with normalized weights w_j,

$$
P_{0,M}^{\mathrm{str}}(a)
=\sum_j w_j\,\delta[a-a_{0,j}^{\mathrm{str}}].
$$

A histogram or KDE may be used only as a numerical representation of this empirical measure.

## 3. Construction from residual microstrain

Choose a reference spacing a_ref appropriate to the reduced normal coordinate and define the structural residual microstrain

$$
\epsilon_0^{\mathrm{str}}(x)
=\frac{a_0^{\mathrm{str}}(x)-a_{\mathrm{ref}}}{a_{\mathrm{ref}}}.
$$

Then

$$
a_0^{\mathrm{str}}(x)
=a_{\mathrm{ref}}[1+\epsilon_0^{\mathrm{str}}(x)].
$$

Equivalently, with lambda_0=a_0^{str}/a_ref,

$$
\lambda_0(x)=1+\epsilon_0^{\mathrm{str}}(x),
$$

and

$$
P_{0,\lambda}^{\mathrm{str}}(\lambda)
=\int_{\Omega}w(x)\,\delta[\lambda-1-\epsilon_0^{\mathrm{str}}(x)]\,dx.
$$

This is the cleanest P0 construction when a residual-strain map is available from diffraction or another strain-mapping technique.

## 4. Construction from a residual normal-stress field

If a structural residual normal-stress field sigma_0^{res}(x) is available instead, define

$$
q_0^{\mathrm{res}}(x)=\frac{\sigma_0^{\mathrm{res}}(x)}{E}.
$$

For the same local generalized-LJ constitutive relation used by the candidate propagator, determine the stable local equilibrium branch from

$$
\phi'[\lambda_0(x)]=q_0^{\mathrm{res}}(x),
$$

subject to

$$
\phi''[\lambda_0(x)]>0.
$$

Then push the resulting lambda_0(x) field forward spatially:

$$
P_{0,\lambda}^{\mathrm{str}}(\lambda)
=\int_{\Omega}w(x)\,\delta[\lambda-\lambda_0(x)]\,dx.
$$

This mapping is a constitutive construction and inherits the validity limits of the generalized-LJ normal calibration.

## 5. Ideal defect-free reference

For an ideal stress-free structural state with no resolved heterogeneity,

$$
\epsilon_0^{\mathrm{str}}(x)=0,
$$

hence

$$
P_0^{\mathrm{str}}(a)=\delta(a-a_{\mathrm{ref}}).
$$

This is not an assumed PDF. It is the direct consequence of a spatially uniform prepared state.

The local-traction candidate then keeps a delta measure a delta if every characteristic receives exactly the same stress history. Therefore an ideal homogeneous conservative specimen does not acquire a probability tail merely because a named distribution was omitted. Any initial spread must come from a resolved physical heterogeneity, and any cycle-to-cycle structural evolution must come from a physical state-evolution mechanism.

## 6. Diffraction route to structural P0

### 6.1 Preferred route: spatial strain mapping

A direct map of local d-spacing or residual strain is preferred because it produces samples of a_0^{str}(x) or epsilon_0^{str}(x) that can be pushed forward without interpreting an entire diffraction line shape as a probability density.

Synchrotron energy-dispersive X-ray diffraction has been demonstrated for residual-strain mapping in materials including Al with strain accuracy better than 1e-4 (Korsunsky et al., J. Synchrotron Rad. 9, 77-81, 2002, DOI 10.1107/S0909049502001905).

Synchrotron studies of Al single crystals also show that as-grown crystals can contain subgrain structure and that residual strain can be measured after deformation; one reported level was about 1e-4 in the studied specimen (Okada et al., Mechanical Engineering Journal, 2020, article 19-00634). These values are evidence that structural heterogeneity exists and is measurable, not a universal P0 width to be inserted into this model.

### 6.2 Secondary route: corrected line-profile inversion

For Bragg angle theta and X-ray wavelength lambda_X,

$$
2a\sin\theta=n\lambda_X.
$$

Thus

$$
a(\theta)=\frac{n\lambda_X}{2\sin\theta},
$$

and

$$
\left|\frac{d\theta}{da}\right|=\frac{\tan\theta}{a}.
$$

If, and only if, a corrected strain-only angular profile p_theta(theta) can be interpreted as an incoherent mixture of local spacings, its push-forward is

$$
P_0^{\mathrm{str}}(a)
=p_\theta[\theta(a)]\frac{\tan\theta(a)}{a}.
$$

If the experimental profile is parameterized by psi=2theta instead, the Jacobian becomes

$$
\left|\frac{d\psi}{da}\right|=\frac{2\tan\theta}{a}.
$$

This inversion is **not automatically valid for a raw diffraction peak**. Instrumental broadening, finite coherent-domain size, peak overlap, mosaic/orientation spread, anisotropic strain broadening, and other specimen effects must first be separated. Diffraction line-profile literature explicitly warns that instrumental broadening must be removed before size/strain information is interpreted (for example, Scardi et al., J. Appl. Cryst. 51, 831-843, 2018).

## 7. Why an instantaneous thermal P0 is a different object

Let xi=a-a_ref. Near the stable reference state,

$$
U(a)\approx U(a_{\mathrm{ref}})+\frac12K_a\xi^2,
$$

with the current calibration

$$
K_a=\left.\frac{d^2U}{da^2}\right|_{a_{\mathrm{ref}}}
=\frac{EA_0}{a_{\mathrm{ref}}}
$$

because phi''(1)=1.

Under the explicit assumptions of a classical canonical harmonic coordinate, the positional marginal is then derived as

$$
P_0^{\mathrm{th}}(\xi\mid T)
\propto
\exp\left[-\frac{K_a\xi^2}{2k_BT}\right],
$$

with

$$
\mathrm{Var}(\xi)=\frac{k_BT}{K_a}.
$$

For a finite harmonic chain with fixed total length, the local variance is modified by the length constraint; for M identical spacing coordinates,

$$
\mathrm{Var}(\xi_i)=\frac{k_BT}{K_a}\left(1-\frac1M\right).
$$

This Gaussian form is not postulated; it follows only from the stated quadratic-Hamiltonian/canonical assumptions.

However, the same canonical preparation also has a nonzero velocity distribution. Therefore in general

$$
F_0^{\mathrm{th}}(a,c)\ne P_0^{\mathrm{th}}(a)\delta(c).
$$

Consequently an instantaneous thermal P0 is **not compatible with the strict P0-only initialization** unless temperature and a justified initial velocity law are supplied as additional preparation data.

Aluminum phonon measurements confirm that the real thermal state contains temperature-dependent phonon populations and spectral changes; neutron measurements have reported the Al phonon DOS from 10 to 775 K (Kresch et al., Phys. Rev. B 77, 024301, 2008, DOI 10.1103/PhysRevB.77.024301). A thermal-diffuse-scattering electron study reported a mean Al atomic displacement of order 12 pm perpendicular to Bragg planes (Herring, Microscopy 62 Suppl. 1, S99-S106, 2013). Such quantities characterize fast thermal motion and should not be silently reinterpreted as a static structural P0.

## 8. Thermal diagnostic red flag for the present normal calibration

Using the retained parameters

- a_ref = 2.8627442948e-10 m,
- E = 69 GPa,
- A0 = 6.0338e-20 m^2,
- m = 12.19,
- n = 6,

one obtains

$$
K_a=\frac{EA_0}{a_{\mathrm{ref}}}\approx14.5431\ \mathrm{N/m}.
$$

A naive classical single-coordinate harmonic calculation at 300 K gives

$$
\sigma_a\approx16.88\ \mathrm{pm},
$$

or

$$
\sigma_\lambda\approx0.05895.
$$

The current normal-curvature threshold is

$$
\lambda_c\approx1.10777154,
$$

corresponding to an offset of about 30.85 pm from a_ref. If the naive harmonic positional Gaussian were directly combined with that threshold, the instantaneous upper-tail mass would be about 3.38 percent at 300 K.

This number is **not a predicted crack probability**. It is a consistency warning: the current instantaneous thermal coordinate, current reduced LJ calibration, and current first-passage crack threshold cannot be combined naively. Relative-displacement correlations, quantum/phonon physics, coarse-graining, and the physical meaning of the threshold must be resolved first.

Therefore the strict P0-only candidate should presently use P0^{str}, not the naive instantaneous thermal P0^{th}.

## 9. What P0 can and cannot solve

A physically constructed P0 solves the initialization problem:

$$
\text{measured/computed prepared structure}\longrightarrow P_0.
$$

It does not by itself solve the laboratory-fatigue time-scale problem. In particular, the local conservative LJ coordinate retains an atomic natural time scale. At laboratory loading frequencies it approaches a quasistatic reversible response unless another physically derived slow or irreversible mechanism is present.

Therefore the candidate local-traction propagator remains quarantined until both of the following are established:

1. a defensible structural P0 construction for the actual specimen;
2. a defensible mechanism that produces the observed cycle-dependent evolution at laboratory fatigue frequencies.

If the second condition cannot be satisfied, return to the exact finite-chain/correlation-hierarchy checkpoint rather than inserting arbitrary diffusion, damping, or a fitted lifetime distribution.

## 10. Current recommendation

For the present research program, define the initialization workflow as

$$
\epsilon_0^{\mathrm{str}}(x)\ \text{or}\ a_0^{\mathrm{str}}(x)
\longrightarrow
P_0^{\mathrm{str}}(a)
\longrightarrow
\mathcal T_{t,0}^{\sigma}[P_0^{\mathrm{str}}].
$$

The first experimental target should therefore be a spatially resolved residual-spacing / residual-strain measurement on the prepared single-crystal Al specimen, preferably with enough resolution to build an empirical push-forward P0 rather than fitting a named PDF.
