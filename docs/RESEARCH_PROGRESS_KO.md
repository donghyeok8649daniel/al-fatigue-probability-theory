# 연구 진전 설명 — 수식 없이 읽는 현재 결론

## 이번에 실제로 새로 확정한 것

주기하중을 받는 흡수형 확률방정식을 한 cycle 동안 그대로 계산하면
“한 주기 생존연산자”가 생긴다. 이 연산자는 임의로 만든 피로식이 아니라
기존 Smoluchowski 방정식의 시간적분 결과다.

이 연산자에는 가장 오래 살아남는 고유분포와 그 분포가 한 cycle 뒤에
얼마나 남는지를 나타내는 최대 고유값이 있다. 이 최대 고유값 하나로
장시간 cycle별 생존비, cycle당 개시확률, 누적 hazard가 동시에 정해진다.
별도 Weibull 수명식이나 damage 계수를 붙이지 않았다.

초기분포가 고유분포와 달라도 두 번째 고유값이 초기조건의 영향이 얼마나
빨리 사라지는지 결정한다. 왼쪽 고유함수는 초기 spacing 위치별 장기 생존
기여도를 정한다. 따라서 장기 생존식 앞에 붙는 초기조건 계수도 fitting하지
않고 계산할 수 있다.

## 왜 기존 예제가 너무 적은 cycle에서 무너졌는가

현재 무차원 예제의 한 cycle 생존비는 약 0.9048이다. 즉 한 cycle마다
intact 확률의 약 9.5퍼센트가 변곡점 흡수경계를 통과한다. 그래서 절반이
개시되는 cycle이 약 6.9회밖에 안 된다.

이는 프로그램이 cycle 수를 잘못 세어서도 아니고 초기분포가 이상해서도
아니다. 두 번째/첫 번째 고유값 비가 약 0.000091이어서 초기분포 영향은
사실상 한 cycle 뒤 사라진다. 낮은 cycle 수는 현재 선택한 무차원 온도,
하중, mobility 시간척도와 개시경계의 직접 결과다. 실제 단결정 알루미늄
수명이라는 뜻은 전혀 아니다.

## 확률분포에서 히스테리시스가 생기는 이유

하중 최고점은 cycle의 25퍼센트 지점인데 평균 spacing, 평균 상호작용
에너지와 hazard의 최고점은 약 35퍼센트 지점, 분산 최고점은 약 40퍼센트
지점에 나타난다.

원인은 finite mobility다. 하중이 바뀌어도 확률분포가 즉시 새 평형분포로
이동하지 못하고 probability current를 통해 유한한 시간 동안 이동한다.
그래서 같은 응력에서도 loading 때와 unloading 때의 분포가 다르다.
diffusion 항 하나가 히스테리시스를 만든다는 설명은 틀리다. drift,
diffusion과 하중시간척도가 함께 정하는 유한시간 수송지연이 원인이다.

## 에너지가 cycle마다 쌓이는가

현재 Markov Smoluchowski 모델에서는 장시간 뒤 생존조건부 분포가 매
cycle 똑같이 반복된다. 따라서 생존한 원자근방의 평균에너지와 tail도
매 cycle 계속 커지지 않는다. 누적되는 것은 변곡점 경계를 넘어 intact
모집단에서 빠져나간 확률이다.

히스테리시스 work는 경로의존 소산량이고 원자간 potential은 상태함수다.
소산 work의 일부를 임의 비율로 potential에 더하면 물리가 아니라 fitted
damage 모델이 된다. 지속적인 조건부 에너지 또는 tail 누적을 원한다면
원자동역학이나 bath 제거에서 유도한 느린 내부좌표 또는 memory kernel이
먼저 있어야 한다.

## 이것이 소성변형인가

normal spacing branch만 보면 소성변형이 아니다. 그러나 이제 별도의 선택적
1D registry branch를 활성화했다. 이 branch는 두 원자열의 정확한
Poisson--Bessel 에너지를 사용하고, registry를 한 주기 안으로 접지 않고
여러 lattice well에 걸쳐 펼친 확률밀도를 직접 계산한다.

소성 판정은 단순히 barrier를 한 번 넘었다는 사실이 아니다. 하중을 제거하고
충분한 relaxation 시간을 준 뒤, well 내부 평균변위는 0으로 돌아왔는데도
평균 well 번호가 0이 아니면 잔류 slip으로 정의한다. 현재 무차원 pulse
예제에서는 하중 제거 뒤 well 내부 변위가 약 0으로 복귀했지만 평균 well
번호는 약 0.4913 남았다. 반대로 대칭 zero-mean 하중 6주기에서는 약
0.0010만 남아 방향성 slip이 거의 상쇄됐다.

따라서 Bessel 모델은 더 이상 비활성 archive만은 아니다. 정확한 에너지와
unwrapped 확률동역학은 활성화됐다. 다만 전위선, 전위 증식, forest
hardening, backstress 및 multiple slip은 아직 없으므로 이를 정량적인
알루미늄 crystal plasticity라고 부르지는 않는다. 또한 normal-chain energy와
two-row registry energy는 다른 기하이므로 서로 더하지 않는다.

## 실제 알루미늄 수명으로 바꾸려면

loading-axis 탄성계수를 사용하면 무차원 normal force는 응력을 그
탄성계수로 나눈 값이 된다. 대표면적은 이 응력비에서는 소거된다.

그러나 thermal 분포 폭에는 대표면적이 남고, 실제 초 단위 시간축에는
spacing friction 또는 mobility가 남는다. 따라서 단결정 방향과 응력만
알아서는 수명을 정할 수 없다. 대표 원자층 면적, 온도, mobility와 변곡점
최초통과 정의를 독립적으로 정하거나 검증해야 한다. FEM element 면적이나
mesh 수를 이 값 대신 쓰면 안 된다.

## 검증한 수치 성질

- 확률밀도는 음수가 되지 않는다.
- 흡수경계에서 잃은 질량과 누적 outgoing flux가 기계정밀도로 같다.
- 작은 격자에서 직접 만든 주기행렬의 최대 고유값과 power iteration이
  일치한다.
- 고유분포에서 직접 여러 cycle을 계산한 생존확률과 고유값의 거듭제곱이
  일치한다.
- cycle 시작 phase를 바꿔도 최대 고유값이 같다.
- grid와 timestep을 줄이면 같은 고유값으로 수렴한다.
- 매우 빠른 cycle에서는 cycle당 유출이 주기에 비례한다.
- 매우 느린 cycle에서는 평균 hazard가 고정하중 subgenerator의 평균
  escape rate로 수렴한다.

## 재현 명령

```powershell
py -3 -m pytest -q
py -3 -m pytest libraries/shear/tests -q
py -3 -m pytest tests/test_registry_plasticity.py -q
$env:MPLBACKEND='Agg'
py -3 -m simulations.run_smoluchowski_floquet
py -3 -m simulations.run_registry_plasticity
```

그래프는 `results/figures/smoluchowski_floquet/`, 원자료는
`results/data/smoluchowski_floquet/`에 있다.
활성 소성 예제는 `results/figures/registry_plasticity/`와
`results/data/registry_plasticity/`에 있다.
