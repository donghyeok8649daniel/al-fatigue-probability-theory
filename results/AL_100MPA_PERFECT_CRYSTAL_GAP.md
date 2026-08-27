# 100 MPa cyclic loading versus the perfect-crystal Al shear scale

This note is a quantitative sanity constraint, not a fatigue-life prediction.

For a uniaxial stress amplitude

$$
\sigma_a=100\ \mathrm{MPa},
$$

the maximum possible Schmid factor for a single slip system is

$$
m_s\le\frac12.
$$

Therefore the resolved shear amplitude is bounded by

$$
\boxed{
\tau_a\le m_s\sigma_a=50\ \mathrm{MPa}.
}
$$

Ogata, Li, and Yip reported a first-principles ideal pure shear strength for Al of approximately

$$
\tau_{\rm ideal}\approx2.76\ \mathrm{GPa}
$$

for the perfect-crystal shear problem considered in their Science 2002 study (DOI: 10.1126/science.1076652).

The scale separation is therefore

$$
\boxed{
\frac{\tau_{\rm ideal}}{\tau_a}
\gtrsim
\frac{2.76\ \mathrm{GPa}}{0.050\ \mathrm{GPa}}
\approx55.2.
}
$$

## Consequence

A model that makes a homogeneous perfect Al slip coordinate cross its ideal atomistic barrier at a 100 MPa axial amplitude by simply reducing the barrier parameter is physically unacceptable unless an independently derived local-amplification mechanism is present.

The missing mechanics must reduce this gap through explicitly modeled structure, for example:

- free surfaces or notches;
- pre-existing defects/non-affine strain fields;
- localized multi-slip interactions;
- persistent slip structures;
- thermally distributed microscopic states;
- crack-like or surface-step stress concentration.

This quantitative gap is useful as a falsification test for future reduced models.

---

# 한국어 번역

축응력 진폭이

$$
\sigma_a=100\ \mathrm{MPa}
$$

이고 단일 slip system의 최대 Schmid factor $m_s=0.5$를 사용하면 resolved shear amplitude는 최대

$$
\boxed{\tau_a=50\ \mathrm{MPa}}
$$

이다.

Ogata, Li, Yip의 first-principles 완전결정 Al 계산에서는 ideal pure shear strength가 약

$$
\tau_{\rm ideal}\approx2.76\ \mathrm{GPa}
$$

규모였다.

따라서 두 응력척도의 비는

$$
\boxed{
\frac{\tau_{\rm ideal}}{\tau_a}\approx55.2
}
$$

정도다.

즉 100 MPa의 축 피로응력에서 완전결정의 균일 slip barrier를 직접 넘게 만들기 위해 barrier 파라미터를 임의로 낮추면 안 된다. 실제 이론은 자유표면, 결함장, multi-slip 상호작용, persistent slip 구조, 유한온도 미시상태, 국부 응력집중 등으로 이 약 55배의 간극이 어떻게 줄어드는지를 역학적으로 설명해야 한다.
