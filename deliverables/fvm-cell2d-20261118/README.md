# FVM 최종 검증 자료

2026년 11월 18일 완료한 2차원 coherent_g4 수치 검증 자료다. 8개 격자에서 총 192구간을 계산해 원래 수치 기준을 전부 통과했다.

- [최종 PDF](CRACK_FVM_2차원_최종검증_20261118.pdf)
- [편집용 Word](CRACK_FVM_2차원_최종검증_20261118.docx)
- [실행 코드 ZIP](CRACK_FVM_cell2d_code_20261120.zip)
- [검증 근거](evidence)
- [코드와 실행 방법](../../research/fvm-cell2d/README.md)

Word는 33개 본문 수식을 네이티브 수식편집기 형식으로 제공한다. PDF 9쪽을 검수했다. 문서는 변수 정의, 지배방정식, 초기조건과 경계조건, 시험결과로 구성했다.

11월 20일 공개용 코드 ZIP은 실행에 필요한 코드·입력, 확정 결과와 구성검사를 묶었다. 계산 코드와 입력의 SHA256은 원 검증본과 같다. 공개 JSON의 개인 PC 경로를 정리한 내역과 원본 SHA256은 SOURCE_PROVENANCE.json에 있다. 원래 로컬 납품본의 목록은 evidence/FINAL_DELIVERY.json에 보존했다. 공개 파일의 실제 해시는 PUBLICATION_MANIFEST.json에서 확인한다.

대표 cell 내부의 동조 운동과 지정 하중 아래 cell 사이 독립성을 가정한다. 실제 특성체적·이동도·물리 시간의 보정은 남아 있다.
