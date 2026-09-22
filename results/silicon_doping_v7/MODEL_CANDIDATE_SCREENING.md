# 도핑 원자 에너지 후보의 후속 검토

2026-09-22의 제한된 공개 원자료 탐색이다. 전체 potential family의 부재를
증명한 검색이 아니며, 아래 모델을 다운로드·실행·채택하지 않았다.

| 확인한 항목 | 실제 의미 | 현재 판정 |
|---|---|---|
| NIST의 `1988_Si(B).tersoff` | Tersoff 순수 Si의 Si(B) 매개변수 버전 | B 도핑 Si potential로 사용하면 안 됨 |
| B–Si modified Tersoff, 저에너지 B bombardment 연구 | 200/500 eV implantation용 다원소 상호작용 | 균열 에너지·전하·저응력 장벽 전이성 검증 전 |
| 원소별 Si/B/P potential | 각 원소 자체의 에너지 | 합쳐서 Si–B/Si–P cross interaction을 만들 근거가 아님 |
| 기존 pure-Si GAP | 순수 Si 환경 학습 | B/P/As/Sb 또는 전하 상태 학습을 자동으로 포함하지 않음 |

첫 항목은 이름이 특히 혼동되기 쉽다. NIST가 순수 silicon의 original
parameterization이라고 명시하며 OpenKIM의 species도 Si다.
출처: [NIST Tersoff Si(B)](https://www.ctcms.nist.gov/potentials/entry/1988--Tersoff-J--Si-b/).

두 번째 항목의 출처는 A. Mari Carmen Perez-Martin, Javier Dominguez-Vazquez,
Jose J. Jimenez-Rodriguez, *A MD study of low energy boron bombardment on silicon*,
Nuclear Instruments and Methods B 164–165 (2000), 431–440,
[DOI 10.1016/S0168-583X(99)01166-0](https://doi.org/10.1016/S0168-583X(99)01166-0).
이번에 접근한 초록/방법 일부는 충돌 문제의 적합성을 다루며,
우리의 균열·캐리어 자유에너지 검증을 대신하지 않는다.

후속 수치 연구의 입력은 최소 다음처럼 구분해야 한다:

1. 순수 Si 중성 기준: 기존 재료 오차를 먼저 평가.
2. 같은 supercell의 B/P 치환 중성 기준: 원소 변화와 전자 변화가 함께 존재.
3. 원소 배치를 고정한 과잉 전자 수 대조: 전자 reservoir와 electrostatic 기준 명시.
4. bulk와 표면/균열 끝의 치환 위치 대조 및 크기·농도 수렴.
5. 각 조건에서 에너지뿐 아니라 force/Hessian·반응 경로를 대조.

균열 부근의 전하 국소화에는 bulk 균일 background 가정의 적용을 다시 검토한다.
보유한 SW/Tersoff의 숫자만 농도에 따라 조정하는 방식은 실행하지 않는다.
