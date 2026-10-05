# Si v6 — 원자 에너지와 시편의 MPa 하중 연결

## 목적

GPa를 MPa로 표시하면 수치는 1,000배가 된다. 재료가 약해지지는 않는다.
이번 연구는 **주어진 균열 형상에서 원격 시편 응력이 원자 균열 끝에 전달되는
과정**을 추가한다. 에너지나 힘에 임의 감소계수를 곱하지 않는다.

공통 구조는 원자 에너지 U → 구속된 상태의 에너지/자유에너지 → 검증된
확률 연산자다. Griffith 관계는 에너지와 단위를 대조하는 기준이며,
그 아래에서 확률을 0으로 만드는 강도 cutoff·소재별 새 피로 법칙이 아니다.

생산 Al/SG/PDE, Si 실행 gate, 이동도와 물리 Hz를 변경하지 않는다.
도핑·산화막·표면 화학을 보정한 결과가 아니다.

## 1. 서로 다른 세 응력 척도

- σ∞: 시편의 원격 명목응력, MPa.
- K_I: 균열 끝의 탄성장 세기를 나타내는 응력확대계수, MPa√m.
- T(δ): 원자면을 균일하게 분리할 때의 traction, GPa.

균열의 반길이는 c로 표기한다. 확률모델의 원자 normal separation a와 구별한다.

\[
K_I=Y\sigma_\infty\sqrt{\pi c}.
\]

Y는 경계조건·형상에 따라 결정한다. 현재 예시는 무한 평판의 중앙 관통균열
2c에 대해 Y=1이다. 유한 폭 시편·표면균열·노치·목이 얇은 시편은 실제
탄성장/J 적분으로 Y 또는 에너지 방출률을 구해야 한다. STL만 읽고 Y=1을
자동 할당하지 않는다. 매끈한 노치의 유한 응력집중을 날카로운 균열 K와
혼동하지 않는다.

## 2. 같은 potential의 이방성 탄성

진행축 x=[11-2]/√6, 법선 y=[111]/√3, 전면 z=[1-10]/√2를 사용한다.
frame의 **행**이 cubic 좌표에서의 국소 기저다. v5에서 같은 potential의
내부 원자를 이완해 얻은 C11/C12/C44를 회전한다.

전면 주기 경계는 plane strain이다. 이 배향은 in-plane/anti-plane이 분리되며
in-plane engineering strain은 (εxx, εyy, 2εxy)다. 분리되지 않는 일반 배향을
2성분 구현으로 조용히 처리하지 않고 거부한다.

\[
\boldsymbol\sigma=D\boldsymbol\epsilon,\quad S=D^{-1},
\]
\[
S_{11}p^4-2S_{16}p^3+(2S_{12}+S_{66})p^2-2S_{26}p+S_{22}=0.
\]

상반평면의 두 root pα를 사용한다. 복소 Airy 계수 dα는
Σdα=1, Σpαdα=0을 만족시킨다. 각 root의 변위 계수는

\[
v_\alpha=(S_{11}p_\alpha^2+S_{12}-S_{16}p_\alpha,
 S_{12}p_\alpha+S_{22}/p_\alpha-S_{26})^T.
\]
\[
u=\frac{2K_I}{\sqrt{2\pi}}\operatorname{Re}
\sum_\alpha v_\alpha d_\alpha\sqrt{x+p_\alpha y}.
\]

같은 식을 미분해 ∇u를 구하고 σ=Dε를 계산한다. 앞쪽 y=0,x>0에서
σyy=K_I/√(2πx), σxy=0이고 뒤쪽 양 균열면은 traction-free다.
등방성의 중근은 불안정한 두 root 차로 나누지 않고 정확한 극한식으로 처리한다.
tip의 연속체 특이 응력을 원자 force로 직접 대입하지 않는다.

Å/GPa 계산에서 1 MPa√m=100 GPa√Å, 1 GPa·Å=0.1 J/m²다.

## 3. 에너지 방출률과 분리 에너지

\[
\mathcal G=H_I K_I^2,
\qquad H_I=\tfrac12\operatorname{Re}
\left[i\sum_\alpha (v_\alpha)_y d_\alpha\right].
\]

GPa 단위의 H와 MPa√m 단위 K를 쓰면
\(\mathcal G[\mathrm{J/m^2}]=1000 H_I[\mathrm{GPa^{-1}}]K_I^2\)다.
등방 plane strain에서는 H_I=(1−ν²)/E다.

별도의 원형 경로 적분으로 검산한다.

\[
J=\oint\left[\tfrac12\sigma:\epsilon\,n_x
 -(\sigma n)\cdot\partial_xu\right]ds.
\]

여러 경로 반경과 각도 해상도에서 J=𝒢를 비교한다.

원자 분리 계산은 하나의 주기 계면을 열어 같은 bulk 배열에서
\(W_{sep}=[U(\delta\to\infty)-U(0)]/A\)를 얻는다. 이 값은 **새 표면 두 개**의
에너지다. 실제 surface relaxation/reconstruction을 계산하지 않은 rigid 기준이며
한 표면 γ와 2γ를 혼동하지 않는다. 3/6/9개의 [111] 주기로 깊이를 비교한다.

셀 높이도 변하므로 traction에는 원자 위치와 셀의 미분이 모두 필요하다.
normal cell height L, 방향 n, upper-half 지시자 I_i에 대해

\[
\frac{dU}{d\delta}=\frac{V}{L}\,n^T\sigma n
 +\sum_i (\nabla_iU\cdot n)
       \left(I_i-\frac{r_i\cdot n}{L}\right).
\]

주기 half-slab의 net force만으로 대체하지 않고 실제 에너지 차분과 비교한다.
\(K_G=\sqrt{W_{sep}/H_I}\)는 지정 에너지 기준의 **Griffith 등식**이다.
원자 lattice-trapping K+·열적 전이율·실측 K_Ic가 아니다.

## 4. 원자 경계와 결정 내부 이완

반경 28/40/56 Å의 원통형 원자 영역, 외곽 8 Å 고정층, 전면 4/8반복으로
계산한다. 고정층의 원자만 K 변위장으로 지정하고 내부 원자의 3성분 변위를
독립적으로 이완한다. 결합을 삭제하지 않는다. 같은 반복 위치로 초기화해도
전체 자유 원자 Hessian을 검사하기 전에는 전면 방향 안정성을 가정하지 않는다.

Diamond는 두 sublattice의 상대 변위 u_int가 있다. 내부 이완 C44를 연속체에
쓰면서 외곽 원자를 affine하게만 고정하면 국소적인 optical 불일치가 남는다.
이를 같은 potential의 원시 힘 미분으로 구한다.

\[
H_{uu}u_{int}+B\epsilon=0,\quad u_{int}=-H_{uu}^{-1}B\epsilon.
\]

6개 cubic engineering strain에 대한 응답을 계산하고 국소 strain을 cubic
축으로 회전해 적용한다. A/B sublattice에 −u_int/2,+u_int/2를 주어 평균
병진을 추가하지 않는다. 이 경계는 선형 내부 이완이며, 균열 끝 위치가 함께
이동하는 flexible boundary나 nonlinear far-field 해법은 아니다.

고정 경계에서 자유 원자가 평형이면 envelope theorem에 따라

\[
\frac{d\,\min_{R_f} U(R_f,R_b(K))}{dK}
 =\sum_{i\in b}\nabla_i U\cdot\frac{dR_{b,i}}{dK}.
\]

경계조건을 바꿔 다시 이완한 ±dK 에너지 차분과 대조한다.
고정 원자의 site energy도 U에 포함한다. vacuum을 포함하는 전체 셀의 virial
응력을 명목 시편 응력으로 보고하지 않는다.

## 5. 확률방정식으로 연결할 때

이 경계조건 아래에서 같은 q의 조건부 U/PMF를 새로 계산해야 한다. v2의 임의
3.5 Å seed-grip 장벽을 이번 K로 단위만 바꿔 대입할 수 없다.

- 변위 제어: 고정 R_b에서의 전체 U/조건부 자유에너지 차이를 사용한다.
- 힘 제어: 외력이 하는 일을 포함한 적절한 potential을 유도한다.
- 이미 전체 원자 경계 일을 포함했다면 −𝒢ΔA를 다시 더해 이중 계산하지 않는다.
- 장벽과 mobility·기억 효과의 검증은 별도이며 원자 질량으로 생산 t0를 정하지 않는다.

Griffith 등식 위에서도 국소 최소가 존재할 수 있다. 그 경우 finite-T 확률은
실제 barrier/PMF와 generator가 정하며, 등식만으로 순간 파손 또는 endurance
limit을 선언하지 않는다. v5의 DFT force 오차와 v4 미완료 연구는 그대로 남는다.

## 6. DFT 자료와의 관계

Cambridge 원자료의 decohesion 33구조는 32/16/24원자 묶음마다 면내 셀이
약 4.5–4.8% 면적 폭으로 달라지고 XC field도 없다. 이것들을 고정 면적의
하나의 rigid opening 경로처럼 미분하거나 끝점 에너지로 재료 Gc를 맞추지 않는다.
이번 새 cleavage는 SW/Tersoff 자체 계산이다. 새 DFT/GAP/MD/fit을 수행한 것이 아니다.

## 출처

- Sih, Paris & Irwin (1965), *On cracks in rectilinearly anisotropic bodies*,
  [DOI 10.1007/BF00186854](https://doi.org/10.1007/BF00186854).
- Dontsova & Ballarini (2017), *Atomistic modeling of the fracture toughness
  of silicon and silicon-silicon interfaces*, Int. J. Fract. 207, 99–122,
  [DOI 10.1007/s10704-017-0224-0](https://doi.org/10.1007/s10704-017-0224-0),
  [저자 원문](https://ballarini.cive.uh.edu/wp-content/uploads/2018/04/106.pdf).
  이방성 원자 경계 접근을 참고했으며 문헌 toughness를 입력값으로 복사하지 않았다.
- Buze & Kermode (2021), *A numerical-continuation-enhanced flexible boundary
  condition scheme applied to Mode I and Mode III fracture*,
  [저자 원문](https://wrap.warwick.ac.uk/149515/1/WRAP-numerical-continuation-enhanced-flexible-boundary-scheme-2021.pdf).
- Bartók et al. (2018), [PRX 8, 041048](https://doi.org/10.1103/PhysRevX.8.041048),
  [공개 원자료](https://doi.org/10.17863/CAM.65004). GAP 학습 자료이며 held-out 검증이 아니다.

실제 계산·종료 상태·검증과 실패 이력은
`results/silicon_specimen_v6/COMPLETED_SUMMARY.md`를 따른다.
