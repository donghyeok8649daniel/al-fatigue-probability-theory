# === 한국어 파일 안내 시작 ===
# - 파일 역할: 활성 1D normal layer-LJ 이론 패키지의 진입점과 공개 모듈 범위를 설명한다.
# - 주요 클래스: 없음 또는 외부 선언만 사용
# - 주요 함수/메서드: 없음 또는 외부 선언만 사용
# - 주의: 이 헤더는 코드 탐색용 설명이며, 물리적 가정/근사 여부는 각 함수 docstring과 docs/의 분류 라벨을 따른다.
# === 한국어 파일 안내 끝 ===
"""Active reduced single-crystal spacing--registry fatigue theory.

The governing intrinsic energy is the multiplicity-free multilayer potential
U0(a,s)=sum_k W(k*a,s). Historical normal-chain and two-row modules remain for
verification, but are not added as the active total energy. The model is a
derived mechanism theory, not a calibrated aluminum-life prediction.
"""
