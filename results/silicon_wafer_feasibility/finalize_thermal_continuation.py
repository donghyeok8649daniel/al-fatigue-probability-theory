"""Create a readable report from completed raw checks without changing them."""
import argparse
from datetime import datetime,timezone
import hashlib,json
from pathlib import Path
import xml.etree.ElementTree as ET

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text(encoding='utf-8'))
def dump(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
def suite(path):
    value=ET.parse(path).getroot()
    if value.tag=='testsuites':value=value.find('testsuite')
    result={k:int(value.attrib[k]) for k in ('tests','failures','errors')}
    result.update(duration_seconds=float(value.attrib['time']),raw_xml_sha256=sha(path))
    return result

def main(a):
    package=a.root/'results/silicon_thermal_continuation';actual=read(package/'atomistic_pilot/summary.json')
    audit=read(package/'analysis/summary.json');protocol=read(package/'atomistic_pilot/protocol.json')
    tests=suite(a.tests_xml)
    if not actual['complete'] or tests['failures'] or tests['errors']:
        raise ValueError('complete calculation and passing related suite required')
    restart=[read(package/f'atomistic_pilot/chain_{i}/restart_check.json') for i in range(2)]
    testrecord=dict(related_pytest=tests,actual_restart_checks=restart,
        new_evaluations_independently_replayed=actual['evaluations_recorded'],
        parent_evaluations_independently_replayed=audit['parent_calls_replayed'],
        all_rng_states_paths_acceptances_and_rejected_repeats_replayed=True,
        first_replay_test_setup=suite(a.setup_error_xml),first_replay_fixture_attempt=suite(a.fixture_error_xml),
        first_setup_cause='default Windows pytest temporary folder denied; corrected with a fresh workspace basetemp',
        first_fixture_cause='synthetic parent fixture omitted completed_proposals; fixture metadata completed, actual run unchanged',
        charts_visually_inspected=a.chart_verified,full_solver_UI_suite_repeated=False,
        kBT_fitted=False,equilibrium_certified=False,material_approved=False,production_changed=False)
    dump(package/'TEST_RESULTS.json',testrecord)
    accept=[r['acceptance_fraction'] for r in audit['chains']]
    result=(
        f'새 원자 에너지·힘 {actual["new_model_calls"]}회, 실제 계산 {actual["elapsed_seconds"]:.3f}초를 완료했다. '
        f'각 체인은 총 저장 96회가 되었으며 새 구간 수용률은 {100*accept[0]:.4f}%/{100*accept[1]:.4f}%다. '
        '기록을 연속으로 늘렸지만 같은 대상 분포에 대한 열평형은 여전히 미인증이다. '
        '평균 차이를 실제 자유에너지 차이·균열 장벽·개시 확률로 바꾸어 부르지 않는다.')
    rows=['| 저장 구간 | 출발 구조 | 평균 에너지, eV | 공통 투영, Å | 열항 잔차, eV | 블록 SE 진단, eV |',
          '|---|---|---:|---:|---:|---:|']
    labels={'parent_32':'앞선 32회','extension_64':'추가 64회','last_32':'마지막 32회'}
    for name,label in labels.items():
        for i,r in enumerate(audit['windows'][name]['chains']):
            rows.append(f'| {label} | {"최초 8%" if i==0 else "복귀 8%"} | {r["mean_features"][0]:.6f} | {r["mean_features"][1]:.6f} | {r["mean_identity_residual_eV"]:+.6f} | {r["residual_block_SE_diagnostic_eV"]:.6f} |')
    rhats=['| 저장 구간 | 에너지 | 공통 투영 | 중점 RMS |','|---|---:|---:|---:|']
    for name,label in dict(labels,combined_96='전체 96회').items():
        rhats.append('| '+label+' | '+' | '.join(f'{v:.6f}' for v in audit['windows'][name]['classical_split_rhat'])+' |')
    identities=['| 마지막 32회 | 자기 출발점 기준 잔차, eV | 공통 중점 기준 잔차, eV | 자기 기준 x/y/z 잔차, eV |',
                '|---|---:|---:|---|']
    last=audit['windows']['last_32']
    for i,(r,c) in enumerate(zip(last['chains'],last['common_reference_identity'])):
        identities.append(f'| {"최초 8%" if i==0 else "복귀 8%"} | {r["mean_identity_residual_eV"]:+.6f} | {c["mean_identity_residual_eV"]:+.6f} | '+', '.join(f'{v:+.6f}' for v in r['Cartesian_xyz_mean_residual_eV'])+' |')
    maxforce=max(r['force_error_eV_A'] for r in restart)
    maxscore=max(r['maximum_errors']['score'] for r in audit['chains'])
    validation=(f'관련 시험 {tests["tests"]}개가 {tests["duration_seconds"]:.3f}초에 통과했다. '
        f'재개 지점의 실제 에너지 차이는 {max(r["energy_error_eV"] for r in restart):.3g} eV, 힘의 최대 차이는 {maxforce:.3g} eV/Å다. '
        f'앞선 {audit["parent_calls_replayed"]}개와 새 {actual["evaluations_recorded"]}개 평가의 좌표·힘 변환·난수·경로·수용·거절 반복을 독립 재생했다. '
        f'새 힘 기울기 재계산의 최대 차이는 {maxscore:.3g}다. 검사에 추가 모델 호출은 없었다.\n\n'
        '독립 재생 시험의 첫 실행은 Windows 기본 임시 폴더 접근이 거부되어 준비 오류가 났다. '
        '쓰기 가능한 새 작업 폴더를 지정했고, 다음 실행에서는 합성 부모 시험 자료에 완료 제안 수가 빠진 것을 발견했다. '
        '합성 자료의 메타데이터를 보완한 뒤 통과했다. 실제 원자 실행과 이전 자료는 이 수정으로 바꾸지 않았다.')
    conclusion=(f'마지막 32회에서도 두 출발점의 평균 에너지 차이는 {abs(last["chains"][0]["mean_features"][0]-last["chains"][1]["mean_features"][0]):.6f} eV, '
        f'공통 투영 차이는 {abs(last["chains"][0]["mean_features"][1]-last["chains"][1]["mean_features"][1]):.6f} Å다. '
        '두 체인의 혼합을 인증하지 않는다. 계산 벽에서의 거절이나 표본 저장 오류가 이 차이를 만든다고 결론 내릴 근거도 없다. '
        '자유에너지 장벽과 그 전이율을 확인하려면 별도의 경로·평형 표본 검증이 필요하다. '
        '다음에는 두 이력이 접근하는 상태 영역을 확인하고, 공통 좌표를 사용하는 샘플링 대조와 영역 크기 검사를 사전에 정해 비교해야 한다.')
    text=a.template.read_text(encoding='utf-8')
    for key,value in [('RESULT',result),('WINDOW_TABLE','\n'.join(rows)),('RHAT_TABLE','\n'.join(rhats)),
                      ('IDENTITY_TABLE','\n'.join(identities)),('VALIDATION',validation),('CONCLUSION',conclusion)]:
        text=text.replace('@@'+key+'@@',value)
    if '@@' in text:raise ValueError('unfilled report field')
    (package/'README.md').write_text(text,encoding='utf-8',newline='\n')
    dump(package/'execution.json',dict(finalized_utc=datetime.now(timezone.utc).isoformat(),
        request='continue active Si research; no renewed duration or automation',
        calculation_started_utc=protocol['started_utc'],calculation_finished_utc=actual['finished_utc'],
        actual_calculation_seconds=actual['elapsed_seconds'],no_extra_agents=True,new_automation_created=False,
        previous_automation_remains_paused=True,initial_si_commit='63bfb3305fb5da1f5acb2ce66651fcc2676bb75a',
        new_MD=0,new_DFT=0,kBT_fitted=False,production_changed=False))
    own=[p for p in package.rglob('*') if p.is_file() and p.name!='package_manifest.json']
    own += [a.root/'solver_v1'/name for name in ('silicon_chain_continuation.py','test_silicon_chain_continuation.py')]
    own += [a.root/'results/silicon_wafer_feasibility'/name for name in (
        'continue_force_thermal.py','audit_thermal_continuation.py','test_thermal_continuation_replay.py',
        'finalize_thermal_continuation.py')]
    references={str((a.root/r).relative_to(a.root)).replace('\\','/'):d for r,d in protocol['inherited_source_sha256'].items()}
    references['results/silicon_force_thermal_comparison/package_manifest.json']=protocol['parent_manifest_sha256']
    dump(package/'package_manifest.json',dict(schema=1,files={str(p.relative_to(a.root)).replace('\\','/'):sha(p) for p in sorted(own)},
        referenced_inputs=references,model_sha256=protocol['model_sha256'],private_XML_not_committed=True,
        manifest_excludes_itself=True))
    print(json.dumps(dict(files=len(own)+1,references=len(references),tests=tests,actual_calls=actual['new_model_calls']),indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--template',type=Path,required=True)
    p.add_argument('--tests-xml',type=Path,required=True);p.add_argument('--setup-error-xml',type=Path,required=True)
    p.add_argument('--fixture-error-xml',type=Path,required=True);p.add_argument('--chart-verified',action='store_true')
    main(p.parse_args())
