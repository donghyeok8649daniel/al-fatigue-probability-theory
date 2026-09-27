"""Build a first-reader report from hash-checked evidence, with explicit draft gates.

No model calls occur. The versioned repository contains all text and figure inputs.
The final addendum is written only after examining actual queue outcomes.
"""
from pathlib import Path
from types import SimpleNamespace
from datetime import datetime, timezone, timedelta
from xml.sax.saxutils import escape
import argparse, copy, hashlib, json, re, sys

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO
BASE = REPO / 'results/silicon_initiation_v13/report_inputs/baseline'
R13 = 'results/silicon_initiation_v13/'
sys.path.insert(0, str(REPO))
from results.silicon_wafer_feasibility.build_initiation_report_v11 import Report, W, H, L, TEAL, GRAY, LINE
from reportlab.pdfbase import pdfmetrics


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rich(text):
    text = str(text)
    for old in ('\u2212', '\u2013', '\u2014', '\u2011'):
        text = text.replace(old, '-')
    text = text.replace('\u2248', '약 ')
    translation = str.maketrans('⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺ᐟ', '0123456789-+/')
    return re.sub('[⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺ᐟ]+', lambda m: '<super>' + m[0].translate(translation) + '</super>', escape(text))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--final-addendum', type=Path)
    parser.add_argument('--preview-outcomes', action='store_true', help='Preview frozen outcome pages without claiming a final commit')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--font', type=Path, required=True)
    parser.add_argument('--bold-font', type=Path, required=True)
    parser.add_argument('--as-of', help='Optional KST display timestamp for reproducible layout')
    args = parser.parse_args()
    final = args.final_addendum is not None
    if final and args.preview_outcomes:
        raise ValueError('Final and preview modes are mutually exclusive')
    output = args.output.resolve()
    if not final and 'output' in output.parts:
        raise ValueError('Drafts must stay outside final output/')
    base = json.loads((BASE / 'content.json').read_text(encoding='utf-8'))
    used = {}
    def use(relative):
        path = REPO / relative
        used[relative] = {'bytes': path.stat().st_size, 'sha256': sha(path)}
        return path
    def data(relative):
        return json.loads(use(relative).read_text(encoding='utf-8'))
    for relative, expected in base['sourceFiles'].items():
        if sha(use(relative)) != expected['sha256']:
            raise ValueError('Previously reviewed evidence changed: ' + relative)
    use('results/silicon_wafer_feasibility/build_initiation_report_v11.py')
    sections = copy.deepcopy(base['sections'][:13])
    sections[0]['subtitle'] = '원자 에너지, 도핑·산화 조건과 첫 균열의 확률 법칙 | 2026.09.28'
    sections[0]['paragraphs'][1] += ' 이번 보고서는 같은 8%의 서로 다른 원자 상태에 대한 곡률 검증과 재료 비교의 표본 범위를 추가해 현재 근거를 설명한다.'
    sections[6]['equations'] = ['∂P/∂t = -∇q·J', 'J = -M(q) [P ∇qF + kBT ∇qP]']
    sections[6]['definitions'].append('q와 F는 동일한 좌표 측도로 정의한다. M(q)는 양의 기계 이동도 행렬이다.')
    sections[12]['takeaway'] = '100 MPa 상태의 국소 곡률 검증은 이전 완료 자료다. 이번 계산 횟수에 다시 합산하지 않는다.'
    sections[12]['paragraphs'][1] = sections[12]['paragraphs'][1].replace('다른 고변형 상태의 전체 곡률 검사는 미완료다.', '새 8%·10% 상태의 결과는 다음 쪽에서 따로 제시한다.')
    for i, section in enumerate(sections):
        if section['kind'] == 'chart':
            section['image_path'] = str(BASE / f'figure-{i + 1:02d}.png')
    # Never turn an existing summary into a certificate without its completion flags.
    rows, sources = [], []
    for name, label in [('loading8', '처음 도달한 8%'), ('return8', '10% 후 복귀한 8%'), ('loading10', '이완한 10%')]:
        relative = R13 + 'dense_' + name + '/summary.json'
        if (REPO / relative).exists():
            item = data(relative)
            complete = item.get('complete') is True and item.get('independent_checked_directions') == 6
            rows.append([label, str(item.get('completed_rows', 0)) + '/' + str(item.get('dimension', 648)),
                         f"{item['minimum_curvature_eV_A2']:.8f}" if complete else '미확정',
                         str(item.get('negative_eigenvalues')) if complete else '미확정',
                         '행렬·6방향 완료' if complete else '부분 결과'])
            sources.append(relative)
        else:
            rows.append([label, '미완료', '미확정', '미확정', '최종 결과 없음'])
    sections.append(dict(title='같은 변위의 서로 다른 국소 상태', kind='table',
        table=dict(head=['상태', '검사한 자유도', '최저 곡률', '음의 고유값', '검사 상태'], rows=rows, widths=[1.25, 1, 1.15, .8, 1.3]),
        takeaway='양의 곡률은 해당 상태 주변의 작은 변형을 검사한다. 전이 장벽이나 첫 균열 여부를 결정하지 않는다.',
        paragraphs=[
            '곡률 단위는 eV/Å²다. 두 8% 상태 모두 고정부를 고정한 648차원 Cartesian 자유도에서 검사했다. 처음 8%의 기존 539행은 입력·모델·기저·기울기 해시를 대조하고 세 행을 다시 계산한 뒤, 남은 109행을 추가했다. 복귀 8%는 648행을 새로 계산했다.',
            '처음·복귀 8%의 잔여 기울기 최대치는 각각 약 0.000722와 0.000704 eV/Å다. 정확한 정지점으로 가정하지 않는다. 가장 작은 고유값의 크기가 달라도 두 상태의 고유벡터가 같은 방향이라는 보장은 없으며, 이 비율을 벌크 영률 변화라고 해석하지 않는다.'
        ], sources=sources))
    if (REPO / R13 / 'grip_augmented_replay.json').exists():
        grip = data(R13 + 'grip_augmented/summary.json')
        checked = data(R13 + 'grip_augmented_replay.json')
        if not grip['complete'] or not checked['complete'] or checked['aggregate_sha256'] != sha(REPO / R13 / 'grip_augmented/summary.json'):
            raise ValueError('Grip comparison is not bound to its independent replay')
        names = {'loading8':'최초 8%', 'return8':'복귀 8%', 'loading10':'10%'}
        sections.append(dict(title='길이를 고정할 때와 힘을 유지할 때', kind='table',
            table=dict(head=['상태', '고정 길이 곡률', '일정 축력 곡률', '맞춘 응력 (GPa)'],
                rows=[[names[x['state']], f"{x['fixed_grip_minimum_curvature_eV_A2']:.8f}",
                       f"{x['force_ensemble_minimum_curvature_eV_A2']:.8f}", f"{x['matched_nominal_stress_GPa']:.6f}"] for x in grip['states']],
                widths=[1,1.35,1.35,1.4]),
            takeaway='길이를 자유롭게 하는 축 방향 좌표를 추가해도 검사한 곡률은 양수다. 각 행은 서로 다른 외력을 유지하는 검사다.',
            paragraphs=[
                '곡률 단위는 eV/Å²다. 자유 원자 216개의 648개 좌표에 양쪽 강체 고정부가 축 방향으로 대칭 이동하는 좌표 하나를 추가했다. 일정 축력을 해당 상태의 반력에 맞추고 U - TΔL의 649차원 곡률을 검사했다. 측면 이동·회전·임의 면하중 조건을 모두 검사한 것은 아니다.',
                '추가한 열만 새로 미분하고 기존 648행 세 벌을 재사용했다. 새 모델 평가 27회와 자동미분 3행을 수행했으며, 독립 힘 유한차분으로 추가 열을 검사했다. 별도 원자료 재검증에서 고유값 차이는 최대 9.593×10⁻¹⁴ eV/Å²였다.',
                '자유원자 힘 잔차가 남아 있으므로 정확한 평형점이나 전역 안정성의 증명은 아니다. 고정부 변위를 허용한 현재 위치의 국소 검사이며 새 힘제어 이완, 벌크 영률 또는 첫 균열 결과와 구별한다.'
            ], sources=[R13 + 'grip_augmented/summary.json', R13 + 'grip_augmented_replay.json', R13 + 'FIXED_AND_FORCE_ENSEMBLES_V13.md']))
    if (REPO / R13 / 'same_grip_modes/summary.json').exists():
        modes = data(R13 + 'same_grip_modes/summary.json')
        replay = data(R13 + 'same_grip_modes_replay.json')
        if not modes['complete'] or not replay['complete'] or replay['summary_sha256'] != sha(REPO / R13 / 'same_grip_modes/summary.json'):
            raise ValueError('Mode comparison is not bound to its independent replay')
        sections.append(dict(title='변형 방향을 고정해 곡률을 비교하면', kind='chart',
            image_path=str(use(R13 + 'same_grip_modes/same_grip_modes.png')),
            points=[f"최저 방향의 각도 {modes['subspaces'][0]['principal_angles_degrees'][0]:.2f}°",
                    f"제곱 내적 {100*modes['lowest_direction_squared_overlap']:.4f}% (확률 아님)",
                    '양끝 위치 차이와 곡률 방향의 기하 비교'],
            takeaway='곡률 변화는 변형 방향을 고정하여 해석해야 한다. 가장 낮은 값 하나로 전체 시편의 강성을 판단할 수 없다.',
            paragraphs=[
                '처음 8%의 가장 부드러운 방향을 그대로 적용하면, 그 방향의 곡률은 처음 상태의 0.029470에서 복귀 상태의 0.097748 eV/Å²로 커진다. 복귀 상태 자체의 가장 낮은 곡률 0.022822는 약 19.95° 다른 방향에서 얻어진다. 그래프의 내적과 투영 비율은 기하학적 양이며 개시 확률이 아니다.',
                '처음 상태의 기울기와 Hessian을 그대로 두고 복귀 위치까지 외삽하면 총 에너지 차이가 +222.234825 eV로 나오지만, 저장 원자료의 실제 차이는 -9.719223 eV다. 큰 유한 변위를 국소 이차식으로 연결할 수 없다는 검사다. 이 차이 또는 외삽의 높이를 전이 장벽으로 쓰지 않는다.'
            ], sources=[R13 + 'same_grip_modes/summary.json', R13 + 'same_grip_modes_replay.json']))
    chord = data(R13 + 'chord_geometry/summary.json')
    sections.append(dict(title='두 상태를 잇는 경로가 필요한 이유', kind='chart',
        image_path=str(use(R13 + 'chord_geometry/chord_geometry.png')),
        points=[f"같은 고정부의 최대 위치 차이 {chord['same_fixed_grip_max_difference_A']:.1f} Å",
                f"자유 원자의 RMS 위치 차이 {chord['labelled_free_atom_RMS_displacement_A']:.6f} Å",
                f"직선 보간 중 최저 원자 간격 {chord['exact_continuous_chord_minimum_distance_A']:.6f} Å"],
        takeaway='양끝 상태의 에너지 차이와 양의 곡률만으로 그 사이 장벽을 알 수 없다.',
        paragraphs=[
            '10%를 거쳐 복귀한 8%의 에너지는 처음 8%보다 9.719223 eV 낮다. 두 상태의 원자 번호를 그대로 잇는 직선 경로에서는 한 쌍이 양끝의 2.329663/2.358467 Å에서 중간의 1.710907 Å로 압축된다. 경로 좌표는 시간이 아니다.',
            '그림은 저장 좌표의 기하 분석이며 새로운 에너지 경로 계산이 아니다. 초기 경로를 바꾸고 실제 에너지로 이완한 뒤, 안장점의 힘·곡률과 양방향 도착 상태를 검사해야 한다. 도착 상태가 새 균열인지 재구성인지도 별도로 확인해야 한다.'
        ], sources=[R13 + 'chord_geometry/summary.json', R13 + 'CURRENT_MODEL_AND_PATH_GATE_V13.md']))
    if (REPO / R13 / 'return_zero_replay.json').exists():
        zero = data(R13 + 'return_zero_force/summary.json')
        zero_audit = data(R13 + 'return_zero_replay.json')
        if not zero['complete'] or not zero_audit['complete'] or not zero_audit['same_external_load']:
            raise ValueError('Completed zero-force return and independent audit required')
        states = zero_audit['states']
        sections.append(dict(title='힘을 제거해도 남는 원자 상태의 차이', kind='table',
            table=dict(head=['같은 외부 축력 0에서 비교', '초기 무축력 상태', '재배열 후 무축력'],
                rows=[['기준 길이 대비 변위 (Å)']+[f"{x['extension_A']:+.8f}" for x in states],
                      ['공칭 축응력 잔차 (GPa)']+[f"{x['nominal_stress_GPa']:.3e}" for x in states],
                      ['자유원자 최대 힘 (eV/Å)']+[f"{x['free_force_max_eV_A']:.3e}" for x in states],
                      ['시편 전체 에너지 (eV)']+[f"{x['energy_eV']:.6f}" for x in states]],widths=[2,1.5,1.5]),
            takeaway=f"초기 무축력 길이보다 {100*zero_audit['extension_difference_over_original_force_free_gauge']:.5f}% 긴 배치가 남았다. 현재 근거는 정적 이력 의존성이다.",
            paragraphs=[
                '10%를 거쳐 8%로 복귀한 기존 상태에서 외부 축력을 0으로 바꾸고 원자와 고정부 길이를 함께 이완했다. 새 모델 평가 147회로 축력·자유원자 힘의 사전 기준 10⁻⁵ eV/Å를 충족했다. 초기 무축력 상태는 이전 완료 자료다. 같은 하중을 가해도 두 상태의 변위는 같지 않다.',
                '길이 차이는 2.237173 Å, 에너지 차이는 재배열 후가 +8.150923 eV다. 단순 균일 변형·병진을 제거한 자유원자 RMS 차이도 0.97592 Å 남는다. 이는 실제 소성이나 에너지 소산의 인증이 아니다. 유한온도에서 오래 유지한 계산을 수행하지 않았다.',
                '여러 원자 간 거리 기준에서 고정부 사이 연결은 유지됐다. 연결 여부만으로 작은 균열의 존재나 부재를 판정할 수 없다. 힘 수렴과 전체 곡률의 안정성은 구분하며, 새 상태의 곡률 검증은 후속 결과에 기록한다.'
            ], sources=[R13 + 'return_zero_force/summary.json', R13 + 'return_zero_replay.json', R13 + 'ZERO_LOAD_AND_NONLINEAR_SCOPE_V13.md']))
    if (REPO / R13 / 'local_harmonic/summary.json').exists():
        harmonic = data(R13 + 'local_harmonic/summary.json')
        if not harmonic['complete']:
            raise ValueError('Completed local harmonic diagnostic required')
        thermal = next(x for x in harmonic['temperature_diagnostics'] if x['temperature_K'] == 300)
        sections.append(dict(title='국소 곡률에서 자유에너지로 가는 간격', kind='chart',
            image_path=str(use(R13 + 'local_harmonic/local_harmonic_diagnostic.png')),
            points=[f"300 K 가정의 조화 보정 +{thermal['quadratic_entropy_term_eV']:.6f} eV",
                    f"두 상태의 정적 에너지 차이 {thermal['delta_energy_eV']:.6f} eV",
                    '새 유한온도 MD 또는 자유에너지 표본 추출 없음'],
            takeaway='국소 이차 에너지를 적분하는 근사는 계산했다. 실제 Si의 유한온도 자유에너지는 아직 확보하지 못했다.',
            paragraphs=[
                '동일 고정부·좌표 측도에서 두 8% 상태의 고전적 단일 조화 우물을 비교했다. ΔF = ΔU + (kBT/2) log(det H복귀 / det H최초)이며, 300 K를 넣으면 ΔF = -9.456409 eV다. 여기서 H는 각 648차원 Hessian이고 에너지는 시편 전체 기준이다.',
                '300 K에서 자유원자 전체의 RMS 변위는 최초·복귀 각각 0.23509/0.23762 Å다. 오른쪽 그림은 원자별 RMS의 최댓값이며 같은 온도에서 0.34774/0.31937 Å다. 유한 우물 경계, 비조화성, 양자 핵효과와 여러 상태의 수는 빠져 있다. 온도별 선은 가정한 온도에 공식을 적용한 결과이며 MD 실험이 아니다. 후속 비선형 검사의 범위도 확인해야 한다.'
            ], sources=[R13 + 'local_harmonic/summary.json', R13 + 'FIXED_AND_FORCE_ENSEMBLES_V13.md']))
    if (REPO / R13 / 'anharmonic_scope/summary.json').exists():
        nonlinear = data(R13 + 'anharmonic_probes/summary.json')
        checked = data(R13 + 'anharmonic_replay.json')
        scope = data(R13 + 'anharmonic_scope/summary.json')
        if not nonlinear['complete'] or not checked['raw_replay_passed'] or not scope['complete'] or checked['verified_points'] != 72:
            raise ValueError('Completed nonlinear probes and raw replay required')
        if checked['reported_summary_sha256'] != sha(REPO / R13 / 'anharmonic_probes/summary.json'):
            raise ValueError('Nonlinear source changed after independent replay')
        sections.append(dict(title='직선 방향의 작은 오차가 충분하지 않은 이유', kind='chart',
            image_path=str(use(R13 + 'anharmonic_probes/anharmonic_line_probes.png')),
            points=['두 8% 상태 × 각 3방향 × 6크기 × 양·음 방향 = 72점',
                    '그림은 최대 원자변위 0.01~0.20 Å의 점', '새 모델 평가 74회 (출발상태 대조 2회 포함)'],
            takeaway='선택한 방향의 에너지 오차가 작아도 수직 방향의 힘 오차는 남는다. 전체 조화 우물의 검증과는 다르다.',
            paragraphs=[
                '각 상태에서 가장 부드러운 세 방향을 따라 원자를 움직여 실제 에너지·힘을 계산했다. 이완하거나 전이 경로를 찾은 것은 아니다. 300 K 이차식의 한 좌표 표준편차 2배에서도 에너지 오차는 최초·복귀 각각 최대 0.06521/0.01833 kBT였다. 이 온도는 변위 크기를 정한 가정이다.',
                '같은 여섯 점에서 가장 큰 전체 기울기 잔차는 0.07750/0.07386 eV/Å이며, 그중 시험 방향에 수직인 성분은 0.07743/0.07385 eV/Å다. 같은 점의 선형 기울기 변화 노름 0.06399/0.06235 eV/Å보다 크다. 모두 648개 자유좌표 전체의 노름이다.',
                '여러 모드를 함께 움직이거나 다른 방향으로 이완할 때의 결합을 단일 직선 에너지로 제한할 수 없다. 72점의 원자료·양음 대조를 재검증했지만 유한온도 자유에너지, 우물 경계, 안장점 또는 균열 개시가 인증된 것은 아니다.'
            ], sources=[R13 + 'anharmonic_probes/summary.json', R13 + 'anharmonic_replay.json', R13 + 'anharmonic_scope/summary.json']))
    sample = data(R13 + 'mpa_sampling_scope/summary.json')
    sections.append(dict(title='재료 모델 비교에서 표본의 의미', kind='table',
        table=dict(head=['DFT 자료', '전체/선택 상태 수', '전체 힘 RMSE', '선택 힘 RMSE'],
            rows=[[d['dataset'], f"{d['full']['frames']} / {d['selected']['frames']}",
                   f"{d['full']['baseline_component_RMSE_eV_A']:.6f}",
                   f"{d['selected']['baseline_component_RMSE_eV_A']:.6f}"] for d in sample['datasets']], widths=[1,1.55,1.2,1.2]),
        takeaway='같은 모델도 선택한 상태에 따라 평균 오차가 크게 달라진다. 같은 원자료 집합에서 비교해야 한다.',
        paragraphs=[
            '표는 기존 MACE-MP-0b3의 힘 성분 RMSE이며 단위는 eV/Å다. CP2K 자료에서는 선택 표본의 오차가 전체보다 크지만, QE 자료에서는 반대로 작다. QE의 30개 조성 중 선택된 것은 20개다. 선택 규칙과 199개 원본 해시는 후보 모델의 새 결과를 보기 전에 고정했다.',
            '이는 원자료 CSV의 재분석이다. 전체 공개 자료도 실제 웨이퍼 조건의 무작위 표본은 아니다. MACE-MP-0b3와 MPA-0는 MPTrj 학습 자료를 공유하므로 두 모델의 일치를 독립적인 오차 보증으로 삼지 않는다. 신규 비교의 실제 완료 범위는 다음 결과 표에서 구분한다.'
        ], sources=[R13 + 'mpa_sampling_scope/summary.json', R13 + 'MODEL_COMPARISON_SCOPE_V13.md', 'https://github.com/ACEsuit/mace-foundations']))
    if (REPO / R13 / 'mpa_replay/summary.json').exists():
        mpa = data(R13 + 'mpa_comparison/summary.json')
        replay = data(R13 + 'mpa_replay/summary.json')
        if not mpa['complete'] or not replay['complete'] or not replay['source_array_checks_exact'] or mpa['frames'] != 199 or replay['frames'] != 199:
            raise ValueError('Complete selected MPA comparison and raw replay required')
        groups = {x['group']:x for x in mpa['groups']}
        selected = [('CP2K__Si_O_H_vacuum1','CP2K Si/O/H'),
                    ('CP2K__Si_O_vacuum0','CP2K Si/O'),
                    ('CP2K__small_Si','CP2K 작은 Si'),
                    ('QE__SiO2_stoichiometry','QE SiO2'),
                    ('QE__Si_only','QE 순수 Si'),
                    ('QE__mixed_other_stoichiometry','QE 기타 Si/O')]
        sections.append(dict(title='같은 구조에 새 재료 후보를 적용한 결과', kind='table',
            table=dict(head=['선택된 구조군', '구조 수', '기존 MP-0b3', '후보 MPA-0'],
                rows=[[label,str(groups[key]['frames']),f"{groups[key]['baseline_force_RMSE_eV_A']:.6f}",
                       f"{groups[key]['candidate_force_RMSE_eV_A']:.6f}"] for key,label in selected],widths=[1.7,.7,1.2,1.2]),
            takeaway='여러 구조군의 힘 오차는 줄었지만 악화된 항목도 있다. 선택 표본의 개선으로 재료 모델을 교체하지 않는다.',
            paragraphs=[
                '표는 힘 성분의 RMSE(eV/Å)다. CP2K 135개와 QE 64개, 총 199개에 새 MPA 예측을 수행했고 표준 계산기 대조 9개를 더했다. 기존 모델 예측은 재사용했다. 새 DFT·MD·피팅은 없다. 군별 수치는 원자·성분을 합쳐 계산했으며 각 구조 오차의 단순 평균이 아니다.',
                'QE SiO2군의 합산 오차는 줄었지만 Si 원자의 오차만 보면 0.112741에서 0.117281 eV/Å로 커졌다. CP2K Si192의 같은 조성 에너지 쌍 오차도 0.005199에서 0.007131 eV/atom으로 커졌다. 힘·에너지·원소별 결과를 구분해야 한다.',
                '이 쪽은 사전에 선택한 199개의 결과다. QE 전체 1,159개와 같지 않다. 전체 자료 보완의 실행 여부는 후속 결과에서 따로 기록한다. 공개 DFT 자료에 맞는 힘과 실제 웨이퍼의 첫 균열·확률·수명은 별도의 검증 단계다.'
            ], sources=[R13 + 'mpa_comparison/summary.json', R13 + 'mpa_replay/summary.json', R13 + 'MPA_SELECTED_COMPARISON_V13.md']))
    if (REPO / R13 / 'mpa_qe_full_replay/summary.json').exists():
        full = data(R13 + 'mpa_qe_full/summary.json')
        full_audit = data(R13 + 'mpa_qe_full_replay/summary.json')
        if not full['complete'] or not full_audit['complete'] or full['frames'] != 1294 or full_audit['reused_prediction_arrays_exact'] != 199:
            raise ValueError('Complete full-QE corpus and exact reuse audit required')
        names = {'QE__O_only':'O만 포함', 'QE__SiO2_stoichiometry':'SiO2 조성',
                 'QE__Si_only':'순수 Si','QE__mixed_other_stoichiometry':'기타 Si/O'}
        groups = [x for x in full_audit['groups'] if x['group'] in names]
        sections.append(dict(title='전체 자료로 확대하면 보이는 반례', kind='table',
            table=dict(head=['QE 전체 화학군', '구조 수', '기존 RMSE', '후보 RMSE', '개선 / 악화'],
                rows=[[names[x['group']],str(x['frames']),f"{x['baseline_force_RMSE_eV_A']:.6f}",
                       f"{x['candidate_force_RMSE_eV_A']:.6f}",f"{x['frames_improved']} / {x['frames_worsened']}"] for x in groups],
                widths=[1.4,.7,1,1,1]),
            takeaway='SiO2의 합산 힘 오차는 줄었지만 493개 중 297개 구조에서는 악화됐다. 전체 지표와 개별 반례를 함께 평가해야 한다.',
            paragraphs=[
                '힘 성분 RMSE 단위는 eV/Å다. QE 전체 1,159개를 보완했고 기존 199개 예측은 정확히 재사용했다. 새 구조 1,095개와 표준 대조 9개를 계산했다. CP2K는 선택 135개를 유지하여 총 기록은 1,294개다. O군 1개와 순수 Si군 14개는 동률이라 마지막 열에 포함하지 않았다.',
                '제곱오차 합은 일부 큰 오차의 감소에 크게 영향을 받는다. SiO2의 구조별 RMSE 변화 중앙값은 오히려 +0.004433 eV/Å다. 선택 표본 16개에서 악화했던 Si 원자의 합산 오차도 전체에서는 0.193093→0.125939로 개선됐다. 자료 범위와 집계 방식이 결론을 바꿀 수 있다.',
                '에너지에도 반례가 남는다. Si40O16의 26개 구조에서 같은 조성 쌍의 에너지 오차는 0.011200→0.018963 eV/atom으로 커졌다. 자료에 대한 개선은 후속 재료 후보의 근거이며, 실제 첫 균열·도핑별 수명·물리 시간의 검증은 아니다.'
            ], sources=[R13 + 'mpa_qe_full/summary.json', R13 + 'mpa_qe_full_replay/summary.json', R13 + 'MPA_FULL_QE_COMPARISON_V13.md']))
    addendum = None
    if final:
        addendum = json.loads(args.final_addendum.read_text(encoding='utf-8'))
        for key in ('reviewed_actual_outcomes', 'git_status_verified', 'manifest_verified'):
            if addendum.get(key) is not True:
                raise ValueError('Final gate missing: ' + key)
        if not re.fullmatch('[0-9a-f]{40}', addendum.get('revision', '')):
            raise ValueError('Final source revision required')
        if not addendum.get('sections'):
            raise ValueError('Actual outcome pages missing')
        if addendum['sections'] != data(R13 + 'report_inputs/final_sections.json'):
            raise ValueError('Final pages differ from the frozen repository inputs')
        sections.extend(addendum['sections'])
        revision = addendum['revision']
    else:
        revision = base['revision']
        if args.preview_outcomes:
            sections.extend(data(R13 + 'report_inputs/final_sections.json'))
        else:
            sections.append(dict(title='후속 검증 - 계산 진행 중', kind='text', points=[
            '세 상태의 곡률·축력 검사와 MPA 전체 QE 비교는 완료',
            '큰 산소 포함 구조는 메모리 사전검사에서 중단',
            '무축력 복귀·비선형 검사 완료, 새 무축력 상태의 곡률 계산 중'],
            takeaway='이 쪽은 초안 상태다. 계획을 완료 결과나 재료 채택으로 보고하지 않는다.',
            paragraphs=['계산 큐가 종료한 뒤 실제 완료·실패·부분 결과와 독립 재검증을 읽어 최종 결과 쪽을 구성한다. 공개 DFT에 대한 새 모델 추론, 이전 자료 재분석, 구현 테스트를 따로 집계한다.'], sources=[]))
    sections.extend(copy.deepcopy(base['sections'][13:]))
    sections[-3]['paragraphs'][0] = sections[-3]['paragraphs'][0].replace('관련 구현 시험은 63개 통과했지만', '이전 v12의 관련 구현 시험은 63개 통과했지만')
    sections[-3]['sources'] += [R13 + 'CURRENT_MODEL_AND_PATH_GATE_V13.md']
    sections[-2]['paragraphs'].append('두 8% 상태 사이의 유한 경로가 발견되더라도 먼저 원자 재배열 경로로 분류한다. 첫 균열 B의 검증이 끝나기 전에는 개시 장벽으로 명명하지 않는다.')
    sections[-1]['paragraphs'][0] = sections[-1]['paragraphs'][0].replace('연구 기준은 silicon-wafer-research의 c3a93f5다.',
        ('현재 결과의 소스 revision은 silicon-wafer-research의 ' + revision[:7] + '다.') if final else
        '기존 결과는 c3a93f5 기준이며, 새 v13 결과는 아직 커밋 전 초안이다.')
    stamp = args.as_of or datetime.now(timezone.utc).astimezone(timezone(timedelta(hours=9))).strftime('%Y-%m-%d %H:%M KST')
    mode = '최종 기록' if final else '검수용 초안 · 연구 진행 중'
    extra_files = {}
    class CurrentReport(Report):
        def use(self, path):
            path = Path(path).resolve()
            try:
                relative = path.relative_to(REPO).as_posix()
                use(relative)
            except ValueError:
                extra_files[path.relative_to(ROOT).as_posix()] = {'bytes':path.stat().st_size, 'sha256':sha(path)}
            return path
        def page(self, title, subtitle):
            if self.n:
                self.finish(); self.c.showPage()
            self.n += 1; self.title = title
            self.c.bookmarkPage('p' + str(self.n)); self.c.addOutlineEntry(title, 'p' + str(self.n), 0, False)
            self.c.setFillColor(TEAL); self.c.rect(L, H-41, 33, 3, fill=1, stroke=0)
            self.c.setFont('KoreanBold', 8); self.c.drawString(L+43,H-43,'SILICON / FIRST CRACK INITIATION')
            self.y=H-69; self.p(rich(title), 'title', 8); self.p(rich(subtitle), 'small', 18)
        def finish(self):
            if self.y < 66:
                raise ValueError(f'Overflow page {self.n}: {self.title} y={self.y}')
            self.pages.append(dict(page=self.n,title=self.title,content_bottom=self.y))
            self.c.setStrokeColor(LINE); self.c.line(L,51,W-L,51)
            self.c.setFillColor(GRAY); self.c.setFont('Korean',7)
            self.c.drawString(L,38, f'{stamp} | {mode}')
            self.c.drawRightString(W-L,38,f'{self.n:02d} / {len(sections)}')
        def close(self):
            self.finish(); self.c.save()
    font, bold = args.font.resolve(), args.bold_font.resolve()
    r = CurrentReport(SimpleNamespace(font=font,bold_font=bold,output=output,results=ROOT,as_of=stamp,label=mode,revision=revision))
    r.c.setTitle('실리콘 웨이퍼 첫 균열 개시 모델 - 2026.09.28')
    r.styles['body'].fontSize=10.6; r.styles['body'].leading=17
    for index,s in enumerate(sections,1):
        r.page(s['title'].replace('\n',' '),s.get('subtitle',f'{index-1:02d} | 현재 원리와 검증 범위'))
        if s['kind']=='cover':
            r.box('<b>'+rich(s['takeaway'])+'</b>')
            r.table(['설명 순서','핵심 질문'],[['목적과 법칙','첫 균열을 무엇으로 정의하고 확률을 어떻게 보존하는가'],['재료와 조건','도펀트·전하·산화층이 에너지와 하중에 어떻게 들어가는가'],['현재 결과','어떤 원자 상태와 오차를 실제로 검사했는가'],['사용 범위','실제 개시 확률과 물리 시간까지 무엇이 더 필요한가']],[1,2.7])
        elif s['kind']=='chart':
            r.picture(s['image_path'],225)
            r.p(rich('조건·수치: '+' / '.join(s['points'])),'small',12)
            r.box('<b>'+rich(s['takeaway'])+'</b>')
        elif s['kind']=='formula':
            for eq in s['equations']:
                formatted = rich(eq).replace('∇q', '∇<sub>q</sub>').replace('kBT', 'k<sub>B</sub>T')
                formatted = re.sub(r'\b([A-Za-z]+)_([A-Za-z]+)\b', r'\1<sub>\2</sub>', formatted)
                r.p(formatted,'section',13)
            for definition in s['definitions']:r.p(rich(definition),'body',9)
            r.box('<b>'+rich(s['takeaway'])+'</b>')
        elif s['kind']=='table':
            t=s['table']; r.table([rich(x) for x in t['head']],[[rich(x) for x in row] for row in t['rows']],t['widths'])
            r.box('<b>'+rich(s['takeaway'])+'</b>')
        else:
            for point in s.get('points',[]):r.p(rich(point),'section',11)
            r.box('<b>'+rich(s['takeaway'])+'</b>')
        r.section('설명')
        for para in s['paragraphs']:r.p(rich(para))
        r.section('근거 자료')
        links=[]
        for relative in s['sources']:
            if relative.startswith('https://'):
                url=relative; label=relative
            else:
                use(relative); label='/'.join(relative.split('/')[-2:])
                if not final and relative.startswith(R13):
                    links.append(rich('로컬 연구 자료: '+label)); continue
                url='https://github.com/donghyeok8649daniel/al-fatigue-probability-theory/blob/'+revision+'/'+relative
            links.append(f'<link href="{escape(url)}" color="#087e8b">{rich(label)}</link>')
        r.p('<br/>'.join(links) if links else '최종 결과의 원자료를 확인한 뒤 추가한다.','small',0)
    r.close()
    # Check every rendered text code point against the registered font cmap.
    all_text = json.dumps(sections,ensure_ascii=False)
    clean = re.sub('<[^>]+>', '', rich(all_text))
    cmap=pdfmetrics.getFont('Korean').face.charToGlyph
    missing=sorted({c for c in clean if ord(c)>31 and ord(c) not in cmap})
    if missing: raise ValueError('Unsupported glyphs: '+repr(missing))
    manifest=dict(as_of_KST=stamp,stage='final_pending_visual_review' if final else 'draft_not_for_delivery',
        source_revision=revision,new_sources_committed=final,source_files=used,additional_files=extra_files,pages=r.pages,
        pdf_sha256=sha(output),builder_sha256=sha(Path(__file__)),baseline_content_sha256=sha(BASE/'content.json'),
        final_addendum_sha256=sha(args.final_addendum) if final else None,font_sha256=sha(font),bold_font_sha256=sha(bold),
        visual_QA_required=True,report_generation_new_model_calls=0,missing_font_glyphs=missing)
    output.with_suffix('.manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    portable_sections = copy.deepcopy(sections)
    for section in portable_sections:
        if 'image_path' in section:
            section['image_path'] = Path(section['image_path']).relative_to(REPO).as_posix()
    output.with_suffix('.sections.json').write_text(json.dumps(portable_sections,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(output=str(output),pages=r.n,minimum_content_bottom=min(p['content_bottom'] for p in r.pages),stage=manifest['stage']),ensure_ascii=False))


if __name__=='__main__': main()
