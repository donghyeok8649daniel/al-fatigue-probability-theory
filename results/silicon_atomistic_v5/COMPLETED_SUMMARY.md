# Si v5 — 원자모델 기반 감사 완료

## 결론

사용자의 “원자모델링부터” 지시에 따라 기존 original SW와 별도의 ordinary
Tersoff 1989 기준을 **같은 원자 배열에서 실제로 계산**했다. 원자 에너지·힘의
수치 구현은 일치하지만, Si 표면·균열의 DFT 힘과 큰 차이가 남는다.
Tersoff로 교체하는 것만으로 재료 보정이 끝나지 않는다. 두 모델 모두
실제 Si 균열/피로 예측용으로 채택하지 않는다.

- 실제 새 정적 계산: 공개 DFT 구조 **2,475개 × 2모델 = 4,950건**.
- 원자료 원자 수 합계 **171,815개**, 구조 종류 27개; 독립 통계 표본 수라는 뜻이 아니다.
- 결정 자체 평형·내부 이완·탄성, 기존 재배열 4구조, 잔류 힘과 fixed-q Hessian 검사 완료.
- 새 MD, 새 DFT, GAP 실행, potential 적합, physical-Hz 보정은 수행하지 않았다.
- 생산 Al/LJ/Bessel/SG/PDE/UI와 Si 실행 gate는 유지했다.

## 1. 원자 에너지의 구현·검증

[`silicon_atomistic_reference.py`](../../solver_v1/silicon_atomistic_reference.py)는
source checksum을 요구하는 optional serial LAMMPS adapter다. 원자 위치·셀의
좌표계, 원자 ID, total/site energy, Cartesian gradient, tensile stress 단위를
명시한다. 원자 질량은 run-0에 필요한 입력이며 시간을 보정하지 않는다.

실제 엔진: LAMMPS **22 Jul 2025**, Python/라이브러리 SHA는
[`atomistic_controls.json`](atomistic_controls.json)에 기록했다. 선택적 ASE는
**3.26.0**이다. 시스템 설치/전역 PATH 변경 없이 기존 외부 엔진과 작업 cache의
ASE를 이용했다. 저장소 기본 requirements나 UI 의존성을 늘리지 않았다.

수치 대조 결과:

- 방향 에너지 미분/힘: 최대 절대 오차 2.19×10⁻⁹ eV/Å 이하.
- 회전·병진·셀 복제: 힘 차이 1.22×10⁻¹³ eV/Å 이하.
- 응력과 affine 에너지 미분: 1.81×10⁻⁷ eV 이하.
- SW 독립 direct-triple과 1,152원자 균열 구조의 힘: 4.40×10⁻¹³ eV/Å 이하.
- 실제 DFT 원자료 8구조 × 2모델 독립 대조: 에너지 차이 1.82×10⁻¹² eV,
  힘 차이 1.30×10⁻¹² eV/Å 이하. 원래 왼손 셀도 독립 구현에 그대로 넣어 대조했다.

이 수치는 **동일 에너지의 구현 일치**다. 아래의 DFT 오차와 구분해야 한다.

## 2. 내부 원자 이완이 탄성을 바꾼다

0 K, 각 모델 자체 diamond 평형, engineering shear γ 기준이다.
Strain 간격 10⁻³, 5×10⁻⁴, 2.5×10⁻⁴를 비교했다.

| 모델 | 격자 Å | 에너지 eV/atom | C11 GPa | C12 GPa | C44 affine GPa | C44 내부 이완 GPa |
|---|---:|---:|---:|---:|---:|---:|
| Original SW | 5.430949777 | −4.336600000 | 151.42449 | 76.42212 | 109.75649 | 56.44934 |
| Tersoff 1989 | 5.432114753 | −4.628872714 | 142.51243 | 75.36411 | 118.79442 | 69.01294 |

내부 변위 3개의 힘 방정식을 풀고 3×3 Hessian의 양의 곡률을 확인했다.
마지막 두 strain 간격 간 탄성계수 최대 변화는 0.000178 GPa 미만이다.
초기 에너지 기반 BFGS는 잔류 힘 1.43×10⁻⁷ eV/Å에서 정밀도에 걸렸고,
힘 방정식을 직접 푸는 방식으로 수정한 후 전체 계산을 다시 수행했다.

**주변 원자의 이완을 허용하느냐가 에너지 곡률을 크게 바꾼다.** 이는 q 좌표를
정할 때 affine 원자 에너지·0 K 이완 에너지·finite-T 자유에너지를 구분해야
한다는 직접적인 계산 예다. 실온 실험값과의 fit 성공을 뜻하지 않는다.

## 3. 공개 DFT와의 실제 비교

자료: [Bartók et al., PRX 2018](https://doi.org/10.1103/PhysRevX.8.041048)의
[Cambridge 모델·학습 자료](https://doi.org/10.17863/CAM.65004).
`Si_PRX_GAP.zip`의 `gp_iter6_sparse9k.xml.xyz`를 읽었다. Archive SHA-256:

```text
1d3efe976c53bcd2889c0556f2d600c26bf4603522b5c92fd8633dff7b7ab5c2
```

실제 파일에는 2,475구조가 있다. 별도 데이터 미러의 2,231이라는 설명을
현재 파일의 개수로 쓰지 않았다. GAP **학습 자료**이므로 GAP의 held-out
검증이라고 부르지 않는다. 파일 속 기존 GAP 예측 240개도 이번 계산에 쓰지 않았다.

DFT 방법 표기: **PW91 1,183 / PBE 50 / 미표기 1,242**. PBE 50개는 Pandey
표면 구조다. force 비교는 방법별로 나눴으며, 에너지 offset 기준은 명시적
PW91끼리만 적용했다. 대소문자가 다른 `dft_*` 2,415개와 `DFT_*` 60개도
빠짐없이 읽고 원래 field 이름을 기록했다.

| 구조 | 방법 | 구조 수 | SW 힘 RMSE eV/Å | Tersoff 힘 RMSE eV/Å |
|---|---|---:|---:|---:|
| Diamond | PW91 | 489 | 0.176896 | 0.204127 |
| (111) 균열 | PW91 | 10 | 1.028627 | 0.794703 |
| (110) 균열 | PW91 | 7 | 0.717442 | 0.777878 |
| (111) 표면 | PW91 | 47 | 1.040733 | 0.992666 |
| Pandey (111) 표면 | PBE | 50 | 1.113458 | 1.102912 |
| 분리 구조 | 미표기 | 33 | 0.581828 | 0.533982 |
| Screw dislocation | PW91 | 19 | 0.547574 | 0.617027 |

위 RMSE는 모든 원자의 **Cartesian 성분**을 동등 가중한다. Vector norm RMS와는
√3 차이가 있다. 각 구조를 동등 가중한 값도 JSON에 별도 저장했다. 큐레이션된
관련 구조들이므로 이 수치로 독립 표본의 신뢰구간을 만들지 않는다.

(111) 균열 DFT 힘 자체의 component RMS는 0.608472 eV/Å다. 이 구조 집합에서
두 모델의 오차가 실제 힘 변화 규모보다 크다. 이것은 수치 오차를 넘는 재료
에너지의 차이이며, 장벽·전이율에 곧바로 정량적인 신뢰를 줄 수 없는 이유다.
두 모델의 실패만으로 모든 analytic energy family의 불가능을 주장하지 않는다.

에너지는 frame **548**의 동일 PW91 diamond 구조를 단일 offset 기준으로 썼다.
이 구조는 이완된 완전 결정이 아니다(최대 힘 0.40532 eV/Å). 선택 규칙은 PW91
diamond 프레임 중 최소 E/N이며 potential fit은 없다. 상대 에너지 RMSE:

| PW91 구조 | SW meV/atom | Tersoff meV/atom |
|---|---:|---:|
| Diamond | 28.8054 | 26.1002 |
| (111) 균열 | 141.0278 | 124.6964 |
| (111) 표면 | 219.3716 | 224.4927 |

이 값은 동일 원자 배열의 상대 에너지다. 이완된 표면 에너지 J/m²,
cohesive energy, 혹은 finite-T barrier가 아니다. 미표기/PBE 프레임을 PW91
offset에 붙이지 않았고 해당 항목은 null로 남겼다.

원자료의 왼손 주기 셀 11개는 격자 basis의 세 번째 벡터 부호만 바꿨다.
정수 격자 집합은 같고 원자 Cartesian 위치는 그대로다. 변경 프레임·연산을
전부 기록했으며 위치를 반사하거나 거울 구조로 바꾸지 않았다.

## 4. 기존 SW 재배열 구조를 다른 모델에서 확인

같은 1,152원자 시편, 고정 grip, gap **2.642745461 Å**에서 계산했다.
같은 grip이 각 모델의 같은 응력·에너지 방출률을 뜻하지는 않는다.

| 구조 | 같은 SW 배열의 SW ΔU eV | 같은 배열의 Tersoff ΔU eV | Tersoff 재이완 후 ΔU eV |
|---|---:|---:|---:|
| state00 | 0 | 0 | 0 |
| state04 | +0.116375 | −0.182243 | −0.144309 |
| state01 | −0.316798 | −0.891355 | −0.807959 |
| state02 | −2.213423 | −3.239194 | −3.118592 |

마지막 열은 **Tersoff에서 재이완한 state00** 기준이다. 두 에너지 함수의
절대 영점을 서로 빼지 않았다. 고정 원자와 gap은 바뀌지 않았다.

초기 Tersoff 이완에서 state00/04는 optimizer success에도 잔류 힘 기준을
넘었다(8.59×10⁻⁶ / 2.03×10⁻⁶ eV/Å). 그 실패 기록을 보존하고, 별도 Newton
보정으로 네 상태의 잔류 힘을 5.60×10⁻¹² eV/Å 이하로 낮췄다.

고정-gap Hessian의 최소 고유값은 0.120185–0.121261 eV/Å²로 양수다.
차분 간격 10⁻⁴ / 5×10⁻⁵ Å를 같은 보정 후 상태에서 대조했다. 따라서 여러
**조건부 국소 최소**는 Tersoff에서도 남지만 상대 에너지는 바뀐다.
Full-q 안정성, Tersoff 전이 경로·장벽, 실제 Pandey 구조와의 동일성은 미인증이다.

## 5. 검증·재현·보관

- Si 전체 관련 테스트 + 소재 gate: 73 PASS / 4.36 s.
- 왼손 셀 basis 처리 추가 후 변경 모듈: 8 PASS / 0.19 s.
  중복을 제거한 관련 테스트는 **74개**. 전체 Al solver/UI 회귀 재실행과는 다르다.
- 원본 DFT 2,475개·예측 4,950개·분류별 70개 집계 재계산, RMSE 차이 0.
- 별도 direct-triple / ASE Tersoff 검증 16건 통과.
- 초기 압축 NPZ 검증은 반복 압축 해제로 느려져 해당 검증 프로세스만 중단했다.
  배열을 한 번 읽도록 수정 후 원시 검증 전체 재실행 5.28 s / PASS. 물리 계산은 중단하지 않았다.
- PNG 직접 확인. Source/hash·실제 tests·artifact 목록은 `validation.json`과
  `artifact_manifest.json`에 기록한다. 최종 커밋/원격 일치는 최종 Git 확인 결과를 따른다.

재현에 NumPy/SciPy/Matplotlib/pytest, optional **ASE 3.26.0**과 MANYBODY가
있는 serial LAMMPS가 필요하다. 경로는 실행 환경에서 지정한다. 기존 결과를
보존하려면 `--output` / `--calculation`에 새 경로를 준다.

```text
python -m results.silicon_wafer_feasibility.run_atomistic_controls --output .cache/si-v5-reproduce
python -m results.silicon_wafer_feasibility.polish_atomistic_reconstruction --calculation .cache/si-v5-reproduce
python -m results.silicon_wafer_feasibility.run_dft_material_audit --download --output .cache/si-v5-reproduce
python -m results.silicon_wafer_feasibility.validate_atomistic_audit --calculation .cache/si-v5-reproduce
```

Archive는 약 90.8 MB이며 ignored cache에 둔다. 다운로드 주소·archive/member
hash를 저장해 원본을 다시 받을 수 있다. 원자료의 DFT 힘·에너지와 두 모델의
모든 예측은 `dft_predictions.npz`, 프레임별 진단은 `dft_frame_metrics.csv`다.
보고서에는 개인 절대 경로나 인증 정보가 없다.

최초 control 계산 후 independent periodic-basis helper가 추가되어 adapter의
파일 hash가 바뀌었다. 최초 실행 소스는 `source_snapshots/`에 원래 hash와
일치하도록 보존했다. 실행한 파일과 최종 파일이 같았다고 표시하지 않는다.

## 다음 원자모델 연구

1. **GAP 또는 screened bond-order 실제 실행**으로 이번 표면·균열 기준을 대조한다.
   현재 QUIP 실행기는 없으며 [공식 배포 목록](https://github.com/libAtoms/QUIP#binary-installation-of-quip-and-quippy)에
   native Windows wheel은 없다.
   설치 가능 여부 확인과 공개 모델 확보는 GAP 실행을 뜻하지 않는다.
2. 원자료에 없는 현재 v4의 재배열 최소·saddle을 독립 재료 모델로 계산한다.
   DFT/GAP의 기존 학습 구조에서 맞는 것과 이 새로운 경로가 맞는 것을 구분한다.
3. Glide/shuffle slip과 relaxed surface의 독립 경로·작업을 추가하고,
   필요한 경우에만 최소 에너지 확장의 식별성·held-out 오차를 검토한다.
4. 재료 에너지 검증 후 조건부 분포·기억·물리 시간으로 돌아간다.
   v4의 미완료 56 profile/2 multibasin 체인 및 전면8 경계·dt 이슈는 미해결이다.

유도와 계약: [원자 기반 이론](../../solver_v1/SILICON_ATOMISTIC_FOUNDATION_V5.md).
그림: [PNG](atomistic_audit.png), [PDF](atomistic_audit.pdf).
