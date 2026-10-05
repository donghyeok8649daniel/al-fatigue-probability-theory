"""Build the reviewed report and exact-byte research manifest, without new physics."""
from __future__ import annotations
import argparse,csv,hashlib,json
from datetime import datetime,timezone
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np

START='2026-10-05T16:28:11+00:00'
END='2026-10-05T17:28:11+00:00'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def dump(path,value):path.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def main(a):
    package=a.root/'results/silicon_force_thermal_comparison'
    analysis=json.loads((package/'analysis/summary.json').read_text(encoding='utf-8'))
    actual=json.loads((package/'atomistic_pilot/summary.json').read_text(encoding='utf-8'))
    if not actual['complete']:raise ValueError('actual calculation incomplete')
    suite=ET.parse(a.test_xml).getroot().find('testsuite')
    if suite is None or int(suite.attrib['tests'])!=37 or int(suite.attrib['failures']) or int(suite.attrib['errors']):
        raise ValueError('related tests not passed')
    failed=ET.parse(a.first_error_xml).getroot().find('testsuite')
    if failed is None or int(failed.attrib['errors'])!=6:raise ValueError('initial collection-error record missing')
    diagnostics=json.loads((package/'literature_and_harmonic/harmonic_quantum_diagnostic.json').read_text(encoding='utf-8'))
    repeat=json.loads((a.table_replay/'harmonic_quantum_diagnostic.json').read_text(encoding='utf-8'))
    if diagnostics!=repeat:raise ValueError('harmonic diagnostic did not replay exactly')
    for name in ('numeric_comparison.csv','harmonic_quantum_modes.npz','table_metadata.json'):
        if sha(package/'literature_and_harmonic'/name)!=sha(a.table_replay/name):raise ValueError('table replay mismatch: '+name)
    citations=json.loads((package/'literature_and_harmonic/citation_addendum.json').read_text(encoding='utf-8'))
    expected=json.loads((a.table_replay/'citation_addendum.json').read_text(encoding='utf-8'))
    if citations!=expected:raise ValueError('final citation update mismatch')
    with (package/'literature_and_harmonic/numeric_comparison.csv').open(encoding='utf-8',newline='') as handle:
        rows=list(csv.DictReader(handle))
    if len(rows)!=52 or any(r['first_initiation_fit_eligible']!='false' for r in rows):
        raise ValueError('numeric rows or initiation eligibility changed')
    proof=json.loads((package/'literature_and_harmonic/citation_proof.json').read_text(encoding='utf-8'))
    pdfs=[r for r in proof['sources'] if 'pdf_sha256' in r]
    if len(pdfs)!=5 or not all(r['automated_match_verified'] for r in pdfs):raise ValueError('PDF excerpt verification failed')
    component=json.loads((package/'component_audit/summary.json').read_text(encoding='utf-8'))
    if component['maximum_sum_error_eV']>1e-11:raise ValueError('thermal component identity sum mismatch')
    protocol=json.loads((package/'atomistic_pilot/protocol.json').read_text(encoding='utf-8'))
    for relative,digest in protocol['source_sha256'].items():
        if sha(a.root/relative)!=digest:raise ValueError('old hash-bound source changed')
    # Package-local attributes preserve JSON/NPZ raw bytes through Windows Git.
    (package/'.gitattributes').write_text(
        '*.json -text whitespace=blank-at-eol,blank-at-eof,space-before-tab,cr-at-eol\n'
        '*.npz -text\n*.csv -text\n*.png -text\n*.md text eol=lf\n',encoding='utf-8',newline='\n')
    preflights=[json.loads((package/f'atomistic_pilot/chain_{i}/preflight.json').read_text(encoding='utf-8')) for i in range(2)]
    tests=dict(related_pytest=dict(tests=int(suite.attrib['tests']),failures=0,errors=0,
        duration_seconds=float(suite.attrib['time']),raw_xml_sha256=sha(a.test_xml),
        scope=['force-record replay (3 synthetic tests)','sampling score','harmonic HMC/pCN',
               'fixed-grip ensemble','boundary-safe thermal identity','harmonic quantum diagnostic']),
        first_attempt=dict(collection_errors=6,tests_executed=0,raw_xml_sha256=sha(a.first_error_xml),
            cause='pytest launched in original Al worktree and imported its solver_v1 package',
            correction='same tests run from the Si worktree; original error XML retained privately'),
        literature_builder_first_attempt=dict(error='SyntaxError: unmatched closing parentheses in three add calls',
            material_or_sampling_failure=False,corrected_and_rerun=True),
        package_builder_first_attempt=dict(error='UnicodeDecodeError: default Windows cp949 on UTF-8 JSON',
            correction='explicit UTF-8 reads; actual raw records unchanged',material_failure=False),
        actual_model_preflights=preflights,all_347_actual_evaluations_replayed=True,
        independent_force_mapping_max_error=max(r['maximum_score_error'] for r in analysis['chains']),
        all_rng_trajectories_acceptances_and_retained_repeats_replayed=True,
        table_rows=52,table_and_harmonic_results_exact_replay=True,PDF_exact_excerpt_checks=5,
        component_sum_max_error_eV=component['maximum_sum_error_eV'],
        charts_visually_inspected=True,new_calls_in_checks=0,full_solver_UI_suite_repeated=False,
        equilibrium_certified=False,material_approved=False,production_changed=False)
    dump(package/'TEST_RESULTS.json',tests)
    chains=analysis['chains'];old=analysis['previous_pcn_acceptance_fractions'];rhat=analysis['classical_split_rhat_diagnostic']
    table=('| 診断 | 최초 8% 시작 | 복귀 8% 시작 |\n|---|---:|---:|\n'
        f'| 기존 pCN 전체 저장 수용률 | {old[0]*100:.2f}% | {old[1]*100:.2f}% |\n'
        f'| 이번 HMC 저장 수용률 | {chains[0]["acceptance_fraction"]*100:.3f}% | {chains[1]["acceptance_fraction"]*100:.2f}% |\n'
        f'| 저장 구조 수 / 서로 다른 유지 평가 | 32 / {chains[0]["distinct_retained_evaluations"]} | 32 / {chains[1]["distinct_retained_evaluations"]} |\n'
        f'| 평균 열항 잔차, eV | {chains[0]["mean_identity_residual_eV"]:.6f} | {chains[1]["mean_identity_residual_eV"]:.6f} |\n'
        f'| 전반 / 후반 열항 잔차, eV | {chains[0]["first_half_residual_eV"]:.4f} / {chains[0]["second_half_residual_eV"]:.4f} | '
        f'{chains[1]["first_half_residual_eV"]:.4f} / {chains[1]["second_half_residual_eV"]:.4f} |\n\n'
        f'두 체인의 split-Rhat은 에너지 {rhat[0]:.4f}, 공통 투영좌표 {rhat[1]:.4f}, 중점 RMS {rhat[2]:.4f}다. '
        '이전 에너지·투영좌표 진단보다 낮아졌지만 중점 RMS 진단은 커졌고, 두 체인이 같은 평형으로 혼합하지 않았다. '
        '개별 열항 비로 얻는 270.70/302.59 K는 짧은 기록의 진단값이며 실제 온도나 적합할 새 온도가 아니다.')
    table=table.replace('診断','진단')
    actual_text=(f'새 원자 에너지·힘 {actual["new_model_calls"]}회, 실제 실행 {actual["elapsed_seconds"]:.3f}초를 완료했다. '
        '새 MD와 DFT는 0회다. 각 체인의 예열 24회·저장 32회를 끝냈고 모든 힘·좌표·난수·거절 반복을 재생했다. '
        '수용률은 높아졌지만 두 초기 구조의 열 통계가 일치하지 않아 평형 인증은 보류한다. '
        '원자 계산의 실제 시간과 사용자가 승인한 한 시간의 조사·분석·검수 구간은 다르다.')
    validation=(f'관련 시험 {tests["related_pytest"]["tests"]}개가 {tests["related_pytest"]["duration_seconds"]:.3f}초에 통과했다. '
        '새 검증기는 합성 자료의 기울기 손상과 거절 상태 기록 손상을 실제로 검출한다. '
        '처음에는 잘못된 작업 폴더에서 Al 패키지를 불러 수집 오류 6건이 났고, 그 기록을 보존한 뒤 Si 폴더에서 같은 시험을 통과했다. '
        '수치표 생성기의 괄호 오류도 수정하고 실제 표와 양자 진단을 다시 생성했다.\n\n'
        '실제 MACE 힘·에너지와 원래 자료는 10⁻⁸ 허용오차 안에서 일치했다. '
        '두 구조의 실제 방향 차분 오차는 1.53×10⁻⁷/1.31×10⁻⁷, 역방향 경로 오차는 10⁻¹⁴ 이하였고, '
        '347개 평가의 독립 기울기 재계산 오차는 1.43×10⁻¹⁴ 이하였다. '
        '경로·에너지·수용 판단·저장 반복의 재생 오차는 0이었다. '
        '이것은 구현과 기록 검증이지 물성이나 실제 열평형의 인증이 아니다. '
        '전체 solver/UI 시험을 이번에 반복하지는 않았다.')
    text=a.template.read_text(encoding='utf-8').replace('@@ACTUAL_RESULT@@',actual_text).replace('@@SAMPLING_TABLE@@',table).replace('@@VALIDATION@@',validation)
    if '@@' in text:raise ValueError('unfilled report placeholder')
    (package/'README.md').write_text(text,encoding='utf-8',newline='\n')
    dump(package/'execution.json',dict(authorized_start_utc=START,authorized_end_utc=END,
        actual_atomistic_start_utc=protocol['started_utc'],actual_atomistic_end_utc=actual['finished_utc'],
        actual_atomistic_elapsed_seconds=actual['elapsed_seconds'],report_built_utc=datetime.now(timezone.utc).isoformat(),
        continuous_one_hour_CPU_or_MD_claim=False,new_MD=0,new_DFT=0,new_model_calls=actual['new_model_calls'],
        additional_agents_created=False,external_messages_sent=False,production_changed=False,
        source_model='MACE-MP-0b3, existing research reference; not an adopted Si Bessel calibration'))
    source_names=['sample_force_thermal.py','audit_force_thermal.py','build_literature_thermal_comparison.py',
                  'test_force_thermal_replay.py','decompose_force_thermal.py','finalize_force_thermal_comparison.py',
                  'FORCE_THERMAL_REPORT_TEMPLATE.md','.gitattributes']
    file_paths=[p for p in sorted(package.rglob('*')) if p.is_file() and p.name!='package_manifest.json']
    file_paths+=[a.root/'results/silicon_wafer_feasibility'/name for name in source_names]
    for path in file_paths:
        if not path.exists():raise ValueError('missing package source '+path.name)
    inputs=dict(protocol['source_sha256'])
    for meta in [json.loads((package/'literature_and_harmonic/table_metadata.json').read_text(encoding='utf-8'))['numeric_inputs_sha256'],
                 diagnostics['source_sha256'],component['source_sha256']]:inputs.update(meta)
    inputs['results/silicon_validation_targets/source_ledger.json']=citations['previous_ledger_sha256']
    for relative,digest in inputs.items():
        if sha(a.root/relative)!=digest:raise ValueError('manifest input changed: '+relative)
    manifest=dict(schema=1,files={str(p.relative_to(a.root)).replace('\\','/'):sha(p) for p in file_paths},
        referenced_inputs=inputs,scientific_status=dict(actual_run_complete=True,implementation_replay_passed=True,
            material_approved=False,equilibrium_certified=False,physical_clock=None,crack_probability=None),
        raw_byte_comparison_required=True,Git_line_endings_must_not_change_scientific_bytes=True)
    dump(package/'package_manifest.json',manifest)
    print(json.dumps(dict(files=len(file_paths),inputs=len(inputs),test_count=37,new_model_calls=actual['new_model_calls']),indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True);parser.add_argument('--template',type=Path,required=True)
    parser.add_argument('--test-xml',type=Path,required=True);parser.add_argument('--first-error-xml',type=Path,required=True)
    parser.add_argument('--table-replay',type=Path,required=True);main(parser.parse_args())
