# Auxiliary Note — Aluminum gamma-surface constraints

## Status

**Auxiliary shear-reference document. Not the main project direction.**

The project is now explicitly centered on cyclic normal stress, normal interatomic spacing $a_i(t)$, the state density $P(a,t)$, and normal-opening instability. See `MILESTONE2_NORMAL_DEFORMATION.md`.

This document is retained only because shear/slip energetics may later be useful as a secondary comparison or closure variable. It must not be interpreted as a requirement that the main fatigue model use a gamma-surface.

## Literature constraint

First-principles and atomistic calculations of fcc Al generalized stacking-fault energetics provide physically relevant scales for shear disregistry. Representative sources include:

- G. Lu, N. Kioussis, V. V. Bulatov, E. Kaxiras, Phys. Rev. B **62**, 3099 (2000), DOI: 10.1103/PhysRevB.62.3099.
- C. Brandl, P. M. Derlet, H. Van Swygenhoven, Phys. Rev. B **76**, 054124 (2007), DOI: 10.1103/PhysRevB.76.054124.
- S. Ogata, J. Li, S. Yip, *Ideal Pure Shear Strength of Aluminum and Copper*, Science **298**, 807–811 (2002), DOI: 10.1126/science.1076652.

These values are **EMPIRICAL / ATOMISTIC VALIDATION INPUTS**, not fatigue evolution laws.

## Relation to the current project

The previous auxiliary slip model used

$$
V_\gamma(s)=\frac{\Delta_\gamma}{2}
\left[1-\cos\left(\frac{2\pi s}{b}\right)\right].
$$

That model was useful only as a nonlinear conservative proof-of-principle for cycle-state evolution.

The current mainline energy model instead prioritizes the fixed generalized Lennard-Jones pair law

$$
v(r)=\varepsilon_{\rm LJ}
\left[
\left(\frac{\sigma_{\rm LJ}}{r}\right)^m
-
\left(\frac{\sigma_{\rm LJ}}{r}\right)^n
\right]
$$

and derives normal spacing dynamics from the actual interatomic geometry and pair-distance distributions.

If a future normal-deformation closure derivation proves that tangential disregistry is necessary, this document provides a reference for adding that variable without inventing an arbitrary shear barrier.

## Rule

Do not lower a shear barrier merely to make a desired macroscopic stress cause slip. Do not promote shear/slip to the main fatigue mechanism unless the normal-deformation theory itself shows that it is required.

---

# 한국어 번역 — Al gamma-surface 보조 참고문서

## 상태

**전단 관련 보조 참고문서이며 프로젝트의 메인 방향이 아니다.**

현재 프로젝트는 반복 수직응력, 수직 원자간격 $a_i(t)$, 상태밀도 $P(a,t)$, 수직 opening instability를 중심으로 한다. 메인 이론은 `MILESTONE2_NORMAL_DEFORMATION.md`를 따른다.

이 문서는 향후 shear/slip energetics가 secondary comparison 또는 closure variable로 필요할 가능성 때문에 남긴다. 메인 fatigue model이 gamma-surface를 반드시 사용해야 한다는 뜻은 아니다.

## 문헌 제약조건

FCC Al generalized stacking-fault energetics에 대한 first-principles 및 atomistic 계산은 전단 disregistry의 물리적 energy scale을 제공한다. 대표 참고문헌은 다음과 같다.

- G. Lu, N. Kioussis, V. V. Bulatov, E. Kaxiras, Phys. Rev. B **62**, 3099 (2000), DOI: 10.1103/PhysRevB.62.3099.
- C. Brandl, P. M. Derlet, H. Van Swygenhoven, Phys. Rev. B **76**, 054124 (2007), DOI: 10.1103/PhysRevB.76.054124.
- S. Ogata, J. Li, S. Yip, Science **298**, 807–811 (2002), DOI: 10.1126/science.1076652.

이 값들은 **EMPIRICAL / ATOMISTIC VALIDATION INPUT**이지 fatigue evolution law가 아니다.

## 현재 프로젝트와의 관계

이전 보조 slip model은

$$
V_\gamma(s)=\frac{\Delta_\gamma}{2}
\left[1-\cos\left(\frac{2\pi s}{b}\right)\right]
$$

를 사용했다.

이 모델의 역할은 nonlinear conservative dynamics에서 cycle-state evolution이 가능한지 보여주는 원리증명뿐이다.

현재 메인 energy model은 대신 고정된 generalized Lennard-Jones pair law

$$
v(r)=\varepsilon_{\rm LJ}
\left[
\left(\frac{\sigma_{\rm LJ}}{r}\right)^m
-
\left(\frac{\sigma_{\rm LJ}}{r}\right)^n
\right]
$$

를 우선하고, 실제 interatomic geometry와 pair-distance distribution에서 수직 spacing dynamics를 유도한다.

향후 normal-deformation closure 유도에서 tangential disregistry가 반드시 필요하다고 증명될 경우에만, arbitrary shear barrier를 새로 만드는 대신 이 문서를 참고해 해당 변수를 추가한다.

## 규칙

원하는 macroscopic stress에서 slip이 발생하도록 shear barrier를 임의로 낮추지 않는다. normal-deformation theory 자체가 필요성을 보여주지 않는 한 shear/slip을 주 fatigue mechanism으로 승격하지 않는다.
