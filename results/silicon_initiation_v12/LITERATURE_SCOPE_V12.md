# v12 추가 문헌의 적용 범위 감사

2026-09-25 확인. 논문 제목의initiation/nucleation을 초기 무균열 시편의
최초 개시와 자동으로 동일시하지 않는다. 아래 자료를 새 수명 보정에 사용하지 않았다.

| 1차 자료 | 실제로 확인한 범위 | 이번 연구에서의 사용 제한 |
|---|---|---|
| Ogata 외, MRS OPL737,F8.23,DOI10.1557/PROC-737-F8.23 | 출판사 초록: 이미 균열이 있는 Si와 crack front의 물 분자에 대한QM/MD. 권호2002,온라인 게시2011 | 초기 무균열 시편의 검증 자료로 채택하지 않음 |
| Hao·Hossain, PRB100,014204(2019),DOI10.1103/PhysRevB.100.014204 | 출판사 초록: 비정질silica의 ReaxFF nucleation/propagation 및 원자 수준 응력 분포 | 전문 접근 불가. 초기조건/수치 재현을 확인하지 못했으며 결정Si의 개시 경계나 force-field로 전이하지 않음 |
| Colombi Ciacchi 외,J.Phys.Chem.C(2008),DOI10.1021/jp804078n;저자 원고arXiv0904.2091 | 전체 원고: 산화Si(001)의 물 반응과 인장/압축 비교. 초기 표면 화학 과정 관측이며 정량 반응률을 제공하지 않음 | 시간에 따른 첫 균열 개시율로 쓰지 않음. 이 논문의 습윤 계산에는 도펀트를 넣지 않았음 |
| Y.J.Oh 외,Ab initio study of boron segregation and deactivation at Si/SiO2 interface(2012),DOI10.1016/j.mee.2011.04.036 | 출판사 초록/방법 발췌: B 치환·침입·Si 자기침입 결함 및 전하에 따른 계면 편석을 다룸. LDA/VASP,400eV | 명목B 농도만으로 국소 원자배치/활성전하를 확정할 수 없다는 근거. P/As/Sb 장벽이나 균열 강도로 재사용하지 않음 |

출처:

- https://www.cambridge.org/core/journals/mrs-online-proceedings-library-archive/article/abs/moisture-effects-of-crack-initiation-in-nanocrystalline-silicon-a-hybrid-densityfunctionalmoleculardynamics-study/DB6790EC0EE02FBA752C757FC0C83732
- https://journals.aps.org/prb/abstract/10.1103/PhysRevB.100.014204
- https://arxiv.org/abs/0904.2091
- https://publica.fraunhofer.de/entities/publication/ef33a2ab-bac5-4180-ab75-e2944a19fe37
- https://www.sciencedirect.com/science/article/abs/pii/S0167931711004606

이번 선택의 의미: 화학적으로 약해진 최초 표면 사건을 연구하려면 Si-O-H의
국소 환경과 전하 상태를 검증해야 한다. 그러나 현재 공개 단일점 오차 감사와
기하학적 재배열 진단만으로 개시율을 결정하지 않는다. 문헌의 crack-front 반응은
별도 배경 지식으로 남기고 초기 균열 없는 개시 연구의 검증 자료와 구분한다.
