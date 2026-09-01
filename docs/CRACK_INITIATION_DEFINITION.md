# Operational definition of crack initiation

## Plain-language decision

This project stops at crack initiation. It does not presently calculate the
opening or propagation of an already formed crack.

Before initiation, one microscopic neighbourhood is assumed small enough that
its nearby layer spacings are approximately equal. Its representative spacing
is the random coordinate used by the probability density. The exact
Riemann-zeta homogeneous-lattice energy is therefore the active energy of that
local neighbourhood. The isolated half-chain gap energy is not used by the
active pre-initiation solver.

Crack initiation is defined as the first time this local representative
spacing reaches the point where tangent stiffness becomes zero. For the current
exponents this critical stretch is about 1.1077715386.

This is an operational definition: it identifies the instant when the intact,
locally homogeneous description loses stability. It does not claim that a
macroscopic free surface has already formed.

## Probability consequences

Two computations must remain distinct:

1. A reflecting computation may have probability above the critical stretch.
   That mass is an instantaneous instability-tail diagnostic and can return.
2. The active initiation computation places an absorbing boundary exactly at
   the critical stretch. Probability reaching it is removed from the intact
   ensemble and becomes cumulative initiation probability.

The absorbing calculation therefore reports survival, outgoing flux, hazard,
and cumulative initiation. Its intact density contains no mass above the
critical stretch by construction.

The force-dependent effective-potential barrier remains a useful stricter
comparison. It lies beyond the tangent-instability point below critical force
and merges with it only at critical force. It is not the active initiation
boundary.

For periodic loading, the absorbing evolution over one cycle has principal
survival multiplier $r$. In its periodic conditional mode,
$S_N=r^N$, cycle escape is $1-r$, and integrated cycle hazard is $-\log r$.
This is the long-cycle consequence of first passage, not a separate empirical
definition of damage. Survivor-conditioned energy is periodic; it is not
artificially accumulated to create initiation.

## FEM and UI interpretation

Each FEM element supplies only its scalar loading-axis normal-stress history.
Elements with identical histories share one probability solution; element
count is not treated as the number of statistically independent atomic
samples. If an `initiation_elements.csv` channel exists, the UI can display
initiation, survival, or hazard. If it does not exist, the UI refuses to invent
values and asks for declared/calibrated probability parameters.

## 한국어 결정문

현재 연구는 균열이 이미 열린 뒤의 성장해석이 아니라 균열개시 연구다.
균열 직전의 작은 원자 근방에서는 인접한 층간거리가 거의 같다고 보고,
그 국소 대표간격을 확률변수로 사용한다. 따라서 활성 에너지는 균일격자의
정확한 제타함수 에너지다.

국소 대표간격이 접선강성 0인 임계간격에 처음 도달하는 순간을 균열개시로
정의한다. 임계점을 통과한 확률은 intact 분포로 되돌리지 않는다. 이때
survival은 아직 개시되지 않은 확률, 누적 유출은 initiation probability,
순간 유출률을 survival로 나눈 값은 hazard다.

주기하중에서는 흡수형 시간진화를 한 cycle 적용한 연산자의 최대 생존
고유값을 $r$이라 한다. 주기 조건부분포에서 $S_N=r^N$, cycle 유출은
$1-r$, cycle 누적 hazard는 $-\log r$이다. 이는 별도 damage law가 아니라
위 first-passage 정의의 장기 결과다. 생존조건부 에너지는 주기적이며
균열개시를 만들기 위해 임의로 누적하지 않는다.

이 정의는 임계점에서 거시적 자유표면이 이미 완성됐다는 뜻이 아니다.
균열 전 국소 균질모델이 더 이상 안정하지 않다는 조작적 판정이다.
