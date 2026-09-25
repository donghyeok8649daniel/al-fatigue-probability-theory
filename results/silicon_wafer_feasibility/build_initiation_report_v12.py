"""Korean v12 report built from completed evidence and explicit partial states.

Uses the existing layout helpers without changing the historical v11 report.
The final PDF requires all-page rendering and visual QA outside this builder.
"""
from __future__ import annotations
import argparse,hashlib,json,sys
from pathlib import Path
from xml.sax.saxutils import escape
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from results.silicon_wafer_feasibility.build_initiation_report_v11 import Report,W,H,L,TEAL,GRAY,LINE


class Report12(Report):
    def __init__(self,args):
        super().__init__(args)
        self.c.setTitle('Si 웨이퍼 균열 개시 v12 - 8시간 보완 연구')
    def data(self,path):return super().data('silicon_initiation_v12/'+path)
    def prior(self,path):return super().data('silicon_initiation_v11/'+path)
    def maybe(self,path):
        return self.data(path) if (self.args.results/'silicon_initiation_v12'/path).exists() else None
    def picture(self,path,max_height=245):return super().picture('silicon_initiation_v12/'+path,max_height)
    def page(self,title,subtitle):
        if self.n:self.finish();self.c.showPage()
        self.n+=1;self.title=title
        self.c.bookmarkPage(f'p{self.n}');self.c.addOutlineEntry(title,f'p{self.n}',0,False)
        self.c.setFillColor(TEAL);self.c.rect(L,H-41,33,3,fill=1,stroke=0)
        self.c.setFont('KoreanBold',8);self.c.drawString(L+43,H-43,'SILICON / FIRST CRACK INITIATION / v12')
        self.y=H-68;self.p(title,'title',8);self.p(subtitle,'small',18)
    def finish(self):
        if self.y<65:raise ValueError(f'page {self.n} overflow: y={self.y}, {self.title}')
        self.pages.append(dict(page=self.n,title=self.title,content_bottom=self.y))
        self.c.setStrokeColor(LINE);self.c.line(L,51,W-L,51)
        self.c.setFillColor(GRAY);self.c.setFont('Korean',7)
        self.c.drawString(L,38,f'{self.args.as_of} KST | {self.args.label}')
        self.c.drawRightString(W-L,38,f'{self.n:02d} / 15')
    def close(self):
        self.finish()
        if self.n!=15:raise ValueError('expected15pages')
        self.c.save()
        manifest=dict(as_of_KST=self.args.as_of,label=self.args.label,pages=self.pages,
            output_filename=self.args.output.name,pdf_sha256=hashlib.sha256(self.args.output.read_bytes()).hexdigest(),
            source_files=self.used,source_revision=self.args.revision,visual_QA_required=True,
            builder_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            layout_helper_sha256=hashlib.sha256(Path(__file__).with_name('build_initiation_report_v11.py').read_bytes()).hexdigest(),
            font_sha256=hashlib.sha256(self.args.font.read_bytes()).hexdigest(),
            bold_font_sha256=hashlib.sha256(self.args.bold_font.read_bytes()).hexdigest())
        self.args.output.with_suffix('.manifest.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')


NAMES={'force100':'100MPa 힘제어','loading8':'최초8% 인장','return8':'10% 이후8% 복귀','loading10':'10% 후속 이완'}
def number(value,digits=5):return '미확보' if value is None else f'{value:.{digits}f}'


def main(args):
    r=Report12(args)
    validation=r.data('validation.json')
    interruption=r.data('session_interruption.json')
    dense=r.maybe('dense_independent_audit/summary.json') or r.maybe('dense_first_completed_replay/summary.json')
    dense_rows=[] if dense is None else dense['states']
    augmentation=r.maybe('grip_augmented/summary.json')
    localization=r.maybe('hessian_localization/summary.json')
    zero=r.maybe('zero_force_return_replay/summary.json')
    harmonic=r.maybe('local_harmonic_diagnostic/summary.json')
    anharmonic=r.maybe('local_anharmonic_audit/summary.json')
    continuum=r.data('charge_boundary_continuum/summary.json')
    large=r.maybe('large_oxide_321_1200/summary.json')
    large_replay=r.maybe('large_oxide_replay/summary.json')
    delta_large=r.maybe('large_delta_transfer/summary.json')
    comparison=r.maybe('mpa_stratified_comparison/summary.json')
    comparison_replay=r.maybe('mpa_independent_replay/summary.json')
    radial=r.data('oxide_radial_delta/summary.json');angular=r.data('oxide_angular_delta/summary.json')
    structure=r.data('structural_descriptors/summary.json')
    energy=r.data('total_energy_error_audit/summary.json')
    guard=r.data('charge_scale_guard/reproduction.json')
    regression=r.data('charge_v9_regression/baseline_replay.json')
    boundary=r.data('charge_basin_limit/summary.json');closed=r.data('charge_basin_limit/closed_form_replay.json')
    old_force=r.prior('force_controlled_prism_360_tight/summary.json')
    paper=r.prior('thermal_oxide_source.json')
    source_manifest=r.data('source_manifest.json');artifact_manifest=r.data('artifact_manifest.json')
    ledger=r.maybe('execution_ledger.json')
    complete_hessians=sum(1 for name in NAMES if (args.results/'silicon_initiation_v12'/('dense_'+name)/'summary.json').exists()
        and r.data('dense_'+name+'/summary.json').get('complete'))

    r.page('Si 웨이퍼 균열 개시','v12 | 8시간 계획 연구: 완료 결과와 최대 절전 중단 기록')
    r.box('<b>목표: 초기 균열을 넣지 않은 Si 시편에서 첫 균열이 시작되는 과정.</b><br/>이번에는 국소 곡률, 산화 환경의 재료 오차, 하중 기준, 전하와 흡수 경계의 수학을 보완했다. 실제 개시 장벽·확률·수명은 아직 검증되지 않았다.')
    r.table(['보완 항목','실제 확보한 근거'],[
        ['원자 상태의 곡률',f'전체 제한 Hessian {complete_hessians}/4상태 계산. 독립 고유값과 힘·에너지 차분을 별도로 대조.'],
        ['재료 에너지 적합','거리8계수 및 거리+각도17계수, 총12회 적합. 외부 QE SiO2 오차가 모두 증가하여 두 후보 미채택.'],
        ['미완료 계산','8% 상태는539/648행의 부분 Hessian. 큰 산화 구조·MPA 모델·0축력 복귀는 미실행.'],
        ['확률 구현',f'작은 전이율에서 최대{100*guard["old_maximum_mass_loss"]:.4f}% 누출을 허용하던 검사 수정. 관련{validation["tests_passed"]}PASS.']],[1,2.7])
    r.section('성과의 범위')
    r.p('수치 구현과 실제 재료의 정확성을 구분했다. 양의 Hessian은 해당 상태의 제한된 국소 곡률을 말하며 전역 최소나 균열 개시를 뜻하지 않는다. 새 원자 평가·기존 자료 재분석·합성 행렬 검증을 따로 집계했다.')
    r.p('Al/UI/생산 에너지 및 physical-Hz gate를 보존했다. 새 DFT와 MD는0회다. GPa를 MPa로 임의 축소하는 계수나 원자 질량에서 만든 overdamped clock을 넣지 않았다.','small')
    link=f'https://github.com/donghyeok8649daniel/al-fatigue-probability-theory/tree/{args.revision}'
    r.p(f'계획 기간:2026-09-25 02:40:54-10:40:54KST<br/>Windows 최대 절전:05:08:49-18:02:21KST, 사유 Battery. 연속8시간 실행은 완료하지 못했다.<br/>원자료·코드: silicon-wafer-research / <link href="{link}" color="#087e8b">{escape(args.revision[:12])}</link><br/>작성 시각:{escape(args.as_of)}KST. 길이Å=10<super>-10</super>m.','small')

    r.page('01. 실행 중단을 먼저 공개한다','Windows 전원 이벤트와 실제 계산 파일로 확인했다.')
    r.box('05:08 배터리 사유로 최대 절전 진입 → 18:02 복귀.<br/>요청한 마감10:40은 절전 구간 안에 있었다. 재개 시 시간 제한이 작동했고 후속 원자 계산은 시작되지 않았다.')
    r.table(['구분','확인된 사실'],[
        ['전원 이벤트','Kernel-Power42: Battery. Power-Troubleshooter1의SleepTime/WakeTime으로 공백 확인.'],
        ['100MPa 전체 행렬','649차원 완료,649개 고유값 모두 양수.49회 에너지·힘 평가와 독립 검증 완료.'],
        ['8% 인장 행렬','539/648행 보존. summary의50,464.4초는 절전 포함 경과시간이며 연속CPU 계산시간이 아님.'],
        ['후속 작업','복귀8%/10% Hessian,grip 추가 감사,큰 산화 구조,MPA 대조,0축력 복귀,비선형 변위점은 미실행.'],
        ['최종 상태','관련 계산/대기 프로세스 종료 확인. 새 장시간 원자 계산을 시작하지 않고 결과 정리.']],[1.1,2.8])
    r.section('완료된 연구는 별도로 보존했다')
    r.p('전체100MPa 곡률,7상태 구조 재분석,12회 보존적 보정 적합,실제 확률 누출 반례와 수정,기존 회귀 재실행,전하 경계의 유한·연속 정확해 검증,63개 구현 테스트는 완료했다.')
    r.p('처음 시도한Windows PowerShell 실행 정책 오류와,경고 문구 때문에 pytest 요약을 못 읽은 보조 기록기 오류도 보존했다. 두 구현 문제는 복구했지만 배터리 최대 절전의 실행 공백과는 다른 사건이다.','small')
    r.box('재개하려면 전원 연결을 유지할 수 있는 환경에서 남은 계산을 실행해야 한다. 부분 행렬을 전체 Hessian으로 사용하거나 미실행 표본을 평가 완료로 합산하지 않는다.')

    r.page('02. 무엇이 부족했으며 무엇을 채웠나','v11의 실패와 미완료 기록을 보존하면서 진행했다.')
    r.table(['시작할 때의 부족','이번 보완','여전히 필요한 것'],[
        ['Krylov16/64에서 수렴 mode0개','전체 행렬 직접 계산:649차원 완료,648차원 부분','다른 크기·경계·재료에서의 검증'],
        ['표면/내부 재배열 구분 부족','초기 표면96·내부120 자유원자별 구조 진단','지속되는 새 균열의 basin/경로'],
        ['산화 표면의 큰 힘 오차','조성 제외 적합 완료. 큰 구조·다른 checkpoint 미실행','개시 경로에 대응하는 전자구조 검증'],
        ['변위0과 힘0의 혼동 가능성','100MPa 접선 계산, 기존 힘0 기준 확인','0축력 복귀 미실행. 열적 유지·크기수렴 필요'],
        ['전하 평형과 흡수 경계의 관계','단위 불변 보존 검사와 정확한 경계 극한','실제 전하별 에너지·전이율']],[1.1,1.35,1.2])
    r.section('개시와 성장의 경계')
    r.p('과거 crack seed·K 경계·강체 계면 분리는 성장 또는 정적 분리 연구로 유지한다. 이번 실제 시편에는 내부 절단면이나 미리 만든 균열을 추가하지 않았다. 응력 최대값이나 거리이웃 변화만으로 첫 균열을 선언하지 않는다.')
    r.box('계산이 완료됨 → 구현이 검증됨 → 재료가 맞음 → 개시 상태와 경로가 맞음 → 물리 시간이 보정됨은 서로 다른 단계다. 앞 단계의 통과로 뒤 단계를 자동 승인하지 않았다.')
    if ledger:r.p('실행별 완료·실패·예산 중단은 execution_ledger.json과 WORKING_STATUS.md에 보존했다.','small')

    r.page('03. 전체 곡률을 계산하는 방법','동일한 중성 MACE-MP-0b3 에너지의 힘을 자동미분했다.')
    r.box('허용 변위 r = r0 + B q, B<super>T</super>B=I<br/>H = B<super>T</super>(d²U/dr²)B<br/>힘제어 에너지 U-F·Delta에서 외력 항은 선형이므로 Hessian에는0, 평형 기울기에는-F가 들어간다.')
    r.table(['항목','규약'],[
        ['시편','360Si, 자유216원자, grip144원자. bare 외부 표면, 내부 crack seed0.'],
        ['변위제어','자유원자648좌표. 최초8%는539행 부분; 복귀8%와10%의 Hessian은 미실행.'],
        ['힘제어','동일 자유좌표+정규화한 rigid grip 연장1좌표,총649차원.'],
        ['grip 좌표','양 grip의 ±Delta/2 변위의 Cartesian norm은6. 정규화 좌표와 원래 연장의 변환을 보존.'],
        ['독립 검사','가장 낮은4방향과 무작위2방향. 힘 차분2e-4/1e-4Å, 에너지 차분2e-3/1e-3Å.'],
        ['구현 대조','대칭화 전 비대칭을 기록. NumPy eigh와 별도 SciPy evr 고유값 비교.']],[1,2.8])
    r.p('기하학적 Cartesian 곡률이며 질량 가중 진동수나 생산 확률 모델의 이동도가 아니다. 전체 행렬을 직접 구성해 이전 반복법의 미수렴을 우회했지만, 이전 실패본을 삭제하거나 성공으로 바꾸지 않았다.')
    r.p('정상점의 잔여힘도 함께 보고했다. 양의 Hessian을 계산했다는 사실만으로 잔여힘0이나 정확한 최소점 존재의 수학적 인증을 주장하지 않는다.','small')

    r.page('04. Hessian 결과와 독립 재검증','양의 값, 음의 값, 미완료를 원자료 그대로 구분한다.')
    if dense_rows:
        r.table(['상태','차원','최소 곡률(eV/Å²)','음의 mode','힘 차분 오차(eV/Å²)'],[
            [NAMES[row['state']],row['dimension'],f'{row["minimum_curvature_eV_A2"]:.7f}',row['negative_modes'],f'{row["force_FD_max_error_norm_eV_A2"]:.2e}'] for row in dense_rows],[1.35,.5,1.2,.65,1.2])
        if dense and dense['complete']:r.picture('dense_independent_audit/dense_hessian_checks.png',225)
        else:r.picture('dense_first_completed_replay/dense_hessian_checks.png',225)
        r.p(f'독립 스펙트럼 차이 최대{max(x["independent_spectrum_max_error_eV_A2"] for x in dense_rows):.2e}eV/Å². 대칭화 전 비대칭 최대{max(x["raw_asymmetry_max_eV_A2"] for x in dense_rows):.2e}eV/Å². 검사한 방향 수는 각6개다.','small')
        r.p('힘 차분은 전체 Hessian-vector 곱을 비교했다. 에너지 둘째 차분은 큰 절대 에너지의 차감 때문에 더 작은 간격에서 반올림 오차가 커질 수 있어 두 간격을 함께 남겼다.')
    else:r.box('최종 전체 곡률 재검증이 확보되지 않았다. 저장된 부분 행렬은 보존했으며 안정성을 통과했다고 보고하지 않는다.')
    r.box('이 결과의 대상은 특정 크기·표면·고정부를 가진0K 시편이다. 양의 곡률에서도 다른 basin으로 가는 유한 장벽이 있을 수 있으며, 장벽의 위치·높이·열적 도달률은 별도로 구해야 한다.')

    r.page('05. 고정 길이와 고정 힘은 다르다','100MPa에서는 고정 힘 조건의 전체 허용 좌표를 함께 계산했다.')
    r.box('H = [[H_ff, h], [h<super>T</super>, H_gg]]<br/>자유원자 이완을 포함한 grip 곡률 = H_gg - h<super>T</super>H_ff<super>-1</super>h<br/>원래 연장에 대한 축강성은 이 값에6²을 곱한다.')
    tangent=next((row.get('relaxed_tangent') for row in dense_rows if row['state']=='force100'),None)
    if tangent:r.p(f'100MPa 상태의 유한 시편 접선은{tangent["tangent_using_current_grip_distance_GPa"]:.6f}GPa다. 이전0-100MPa 두 점 secant 약112.163GPa와 비교 가능한 같은 작은 시편의 응답이며, bulk Si 영률을 보정한 값은 아니다.')
    if augmentation:
        r.table(['상태','고정 길이 최소 곡률','고정 힘 최소 곡률','시편 접선(GPa)'],[
            [NAMES[x['state']],f'{x["fixed_grip_minimum_curvature_eV_A2"]:.6f}',f'{x["force_ensemble_minimum_curvature_eV_A2"]:.6f}',number(x.get('finite_prism_tangent_GPa'),4)] for x in augmentation['states']],[1.3,1.1,1.1,1])
        r.p(f'기존 고정단 행렬{augmentation["reused_fixed_grip_rows"]}행을 재사용하고 새 grip 미분{augmentation["new_grip_autograd_rows"]}행만 계산했다. 독립 힘 차분으로 교차 열을 확인했으며 새 에너지·힘 평가{augmentation["new_graph_evaluations"]}회다.','small')
    else:
        force100=r.data('dense_force100/summary.json')
        r.table(['100MPa의 같은 상태','차원','최소 곡률(eV/Å²)'],[
            ['grip 연장을 고정한 부분공간',648,f'{force100["fixed_grip_at_same_state_lowest_eight_eV_A2"][0]:.9f}'],
            ['grip 연장을 허용한 힘제어',649,f'{force100["minimum_curvature_eV_A2"]:.9f}']],[1.8,.5,1.2])
        r.box('위 값은 완료한100MPa 전체 행렬과 그 부분 블록에서 얻었다. 다른 상태의 추가 grip 감사는 미완료다. 고정 길이의 양의 곡률만으로 고정 힘의 안정성을 대체하지 않는다.')
    r.section('어느 힘을 고정했는가')
    r.p('v11에서 외부 명목 응력100MPa에 맞추어 이완한 상태를 사용했다. 이번에는 자유원자와 정규화한 grip 연장 좌표를 함께 미분했다. 선형 외력의 Hessian은0이며, 남은 enthalpy 기울기의 최대 성분은2.9321e-6eV/Å였다. 다른 변위제어 상태의 grip 추가 검사는 미실행이다.')
    r.p('표의 곡률 단위는eV/Å²다. 원래 명목 면적과 현재 grip 간격으로 접선을 환산했으며, 지지 조건이 다른 실험 시편에 그대로 대입하지 않는다.','small')

    r.page('06. 원자 재배열을 더 자세히 구분','처음의 표면·내부 라벨을 유지하고 거리 규약의 영향을 조사했다.')
    r.picture('structural_descriptors/spatial_pair_changes.png',240)
    selected=[x for x in structure['graph_cutoff_results'] if x['cutoff_A']==2.8 and x['state'] in ('loading_8pct','loading_10pct','return_8pct')]
    r.table(['상태','초기 쌍 손실','새 가까운 쌍','내부-내부 손실','독립 graph 경로'],[
        [x['state'],x['lost_from_zero'],x['formed_from_zero'],x['lost_core_core'],x['edge_disjoint_paths']] for x in selected],[1.3,.8,.85,1,1])
    r.p('초기에 표면96개·내부120개로 나눈 두 영역에서 재배열이 나타났다. 위 초기 쌍은 이완된0% 상태의2.8Å 규약을 기준으로 한다. 최대유량으로 구한 graph 경로 수는 단위 용량의 연결성 진단이며 실제 파괴면이나 결합 에너지가 아니다.')
    r.box('이웃이3개뿐이면3차원 affine 변환이 정확히 맞아 D² 잔차가0일 수 있다. rank3이면서 잔차 자유도&gt;0인 원자만 지역 통계에 포함하고2.8/4.2Å 이웃을 함께 검사했다.')
    r.p('2.7-3.7Å로 규약을 바꾸면 경로 수와 쌍 손실이 달라진다. 한 cutoff에서 손실0이라는 사실로 개시가 없었다고 인증하지 않는다. strain8/10%는 최초 grip 길이 기준이며 힘0 길이 기준이 아니다.','small')

    r.page('07. 최소 재료 보정 두 후보는 실패','양의 개선 수치만 골라 채택하지 않았다.')
    r.table(['후보','계수 / rank','조성 제외 RMSE 전후(eV/Å)'],[
        ['거리 kernel',f'8 / {radial["all_fit"]["rank"]}',f'{radial["held_composition_baseline_RMSE_eV_A"]:.6f} → {radial["held_composition_candidate_RMSE_eV_A"]:.6f}'],
        ['거리+각도 kernel',f'17 / {angular["all_fit"]["rank"]}',f'{angular["held_composition_baseline_RMSE_eV_A"]:.6f} → {angular["held_composition_candidate_RMSE_eV_A"]:.6f}']],[1.2,.85,2.1])
    r.p('CP2K487상태를5개 조성으로 나누어 조성 하나 전체를 제외했다. 각 후보의5fold+전체fit, 총12회 실제 선형 적합이다. 힘만 적합했고 에너지 offset이나 배율을 맞추지 않았다. 힘은 같은 에너지의 정확한 미분이며 병진·회전·원자순서·주기 기저 대조를 검사했다.')
    r.box('외부 QE SiO2에서 기준0.158227eV/Å가 거리 보정 후0.171635, 각도 포함 후0.187492로 악화됐다. 두 후보 모두 미채택이다. 각도형은 거리형의 외부 실패를 본 뒤 선택했으므로 QE는 재사용 진단이며 손대지 않은 최종 test가 아니다.')
    r.section('원자당 오차 뒤에 가려지는 전체 에너지 규모')
    selected=[x for x in energy['records'] if x['family']=='angular']
    r.table(['제외 조성','원자수','기준/각도 포함 보정 에너지 차 RMSE(eV)'],[
        [x['composition'],x['atoms_per_frame'],f'{x["baseline_pair_difference_RMSE_eV"]:.4f} / {x["candidate_pair_difference_RMSE_eV"]:.4f}'] for x in selected],[1.3,.65,1.8])
    r.p('같은 조성의 구조쌍 차이를 비교한 값이다. 쌍들은 독립 표본도 전이 경로도 아니므로 이 값을 개시 장벽의 오차라고 단정하지 않는다. 두 후보의 실패가 모든 에너지 함수형의 불가능성 증명은 아니다.','small')

    r.page('08. 실제 재현한 확률 누출 검사 결함','작은 전이율에서도 시간 단위를 바꾸면 같은 확률 법칙이어야 한다.')
    r.box('이전 검사에서는 max(1,전체 전이율) 형태의 절대 척도가 매우 느린 generator의 누출을 가렸다.<br/>L = 1e-15 · [[-1,1],[1,-1.1]], t=1e15를 정확한 이전 Git 소스로 실행해 재현했다.')
    r.table(['관측','실제 결과'],[
        ['이전 propagator의 열 질량',f'{guard["old_column_masses"][0]:.10f} / {guard["old_column_masses"][1]:.10f}'],
        ['허용되었던 최대 질량 손실',f'{100*guard["old_maximum_mass_loss"]:.6f}%'],
        ['수정 뒤 결과',escape(guard['new_rejection'])],
        ['수정 방식','열별 전이율 척도로 보존 검사. 정지성과 상세균형은 지역 확률 flux 척도로 검사.'],
        ['추가하지 않은 보정','확률 clipping, 재정규화, 양의 고유값의 인위적 제거 없음.']],[1.3,2.7])
    r.p('시간척도1e-15/1/1e15, 빠른 상태군 옆의 느린 누출, 정지하지만 비가역적인 순환, equilibrium 정규화,0-generator를 검사했다. underflow나 수치적으로 미해결인 propagator는 보정값을 만들어 반환하지 않는다.')
    r.section('기존 결과도 보존되는지 확인했다')
    r.p(f'기존 v9의 전하·흡수·동결 배치 혼합 감사 전체를 실제 재실행했다. CSV4개와summary1개, 총{regression["byte_equal_files"]}개가 이전과 byte 단위로 같았다. 새 guard는 유효한 기존 합성 계산의 수치를 바꾸지 않았다.')
    r.p(f'최종 관련 테스트:{validation["tests_passed"]}PASS. 검증한 소스의 SHA256과 실제 최종 소스의 결합은 manifest 검증에 포함했다. 구현 검증이 실제 Si 전이율이나 물리 초·Hz 보정은 아니다.','small')

    r.page('09. 전하가 빨라도 개시 경계는 별도','전하별 개시 집합이 다르면 경계를 단순 평균할 수 없다.')
    r.picture('charge_basin_limit/charge_basin_limit.png',225)
    r.box('기존 공통 법칙 L = L_slow + nu K에서 B에 처음 들어간 질량을 보존한다.<br/>공간 위치q=0,1,2와 전하0,1을 둔 합성 예제에서 B는(q2,전하1) 하나다. 전하 평형 점유p는 도펀트 농도가 아니다.')
    r.p('q2의 살아남은 전하0에서 빠져나가는 개시 전이는nu*p다. 따라서nu만 크다고 경계가 빠르게 축약되는 것은 아니다. 공간 운동을 끈 경계 생존의 정확한 답exp(-nu*p*t)를27조건에서 확인했다.')
    r.p('평균 최초 도달시간의 닫힌식:<br/>T=3/r + 3(1-p)(nu²+3nu·r+r²) / [nu·p(nu+r)(nu+3r)]<br/>nu→infinity에서T→3/r=7.5이지만 p=0에서는 닫힌 비흡수 집단이 남는다.','small')
    r.p(f'24개 행렬,120시점. 질량 잔차 최대{boundary["maximum_mass_error"]:.2e}, 독립 Padé/Krylov 차이{boundary["maximum_Pade_vs_Krylov_error"]:.2e}. 평균시간 닫힌식의 최대 차이{closed["maximum_absolute_error_model_time"]:.2e} model time.','small')
    r.p('유한 sink를 평균내는 기존 검사를 무한 흡수 경계에 자동 적용하지 않는다. 이 결과는 고정 공간 격자의 합성 극한이며 Si의 실제 전하 경계·전이율을 보정하지 않았다.','small')

    r.page('10. 빠른 전하와 공간 격자를 함께 검증','고정 격자에서 맞는 근사가 연속 공간의 유한율 오차까지 보장하지 않는다.')
    r.picture('charge_boundary_continuum/joint_charge_spatial_resolution.png',240)
    r.p('두 전하가 같은D로[0,L]를 확산한다. 왼쪽은 둘 다 반사, 오른쪽은 전하0만 반사하고 전하1은 흡수한다. 처음은x=0의 평형 전하 점유(1-p,p)다. 실제 Si 입력을 사용하지 않은 정확해 대조다.')
    r.box('연속 평균시간 = L²/(2D) + (1-p)L coth(L sqrt(nu/D)) / [p sqrt(D nu)]<br/>큰nu의 초과 시간은nu<super>-1/2</super>에 비례한다. 고정 격자의 주도 보정은nu<super>-1</super>이며, 공간 간격은sqrt(D/nu) 경계층을 분해해야 한다.')
    r.p(f'실제 희소 선형계{continuum["synthetic_backward_systems"]}개와 연속 방정식·경계조건{continuum["analytic_boundary_checks"]}조건을 검사했다. 격자 닫힌식 대비 최대 상대 차이{continuum["maximum_discrete_formula_relative_error"]:.2e}, backward 잔차 최대{continuum["maximum_backward_system_absolute_residual"]:.2e}다.','small')
    r.p('p=0.2,nu=10000,D=0.4,L=1에서 연속 평균시간은1.31324555다. N=16이면1.26254564(-3.8607%), N=1024이면1.31305791(-0.0143%)였다. 완전 흡수값1.25를 뺀 초과 시간의 오차도 따로 비교했다.','small')
    r.p('두 극한의 최종 흡수값은 같다. 서로 다른 것은 유한 전하율의 주도 보정이다. 이 평균시간 해를 모든 시간의 단일 상수 경계조건이나 실제 Si 전이율로 해석하지 않는다.','small')

    r.page('11. 도펀트별로 달라져야 할 입력','공통 확률 법칙을 유지하며 화학배치와 전하를 명시한다.')
    r.box('고정 화학배치 D별 생존을 먼저 구한다.<br/>S_mix(t)=sum_D w_D S_D(t)<br/>생존 집단의 배치 비중=w_D S_D(t)/S_mix(t). 평균 에너지·평균률을 한 번 대입하는 것과 일반적으로 다르다.')
    r.table(['구분','실제로 필요한 자료','현재 상태'],[
        ['도펀트 원소 B/P/As/Sb','치환·침입·cluster·계면 편석의 구조와 에너지','이번 Si/O/H 감사로 검증되지 않음'],
        ['명목 농도','국소 배치 분포, 활성 분율, 공간적 편석','작은 셀에 항상1원자를 넣어 대표할 수 없음'],
        ['전하 상태','전자수별 에너지, chemical potential, 경계 전하 조건','중성 MACE는 전자수 제어 모델이 아님'],
        ['개시 경계','각 상태의 intact/전이/지속되는 균열 영역','단순 거리 cutoff나 평균 경계로 대체 불가'],
        ['물리 시간','대응 좌표의 실제 전이·완화·이동도','초·Hz·수명 미보정']],[.9,1.5,1.5])
    r.p('도핑이 항상 더 강하게 또는 더 약하게 만든다는 단일 순위를 제시하지 않는다. 같은 원소와 평균 농도라도 배치, 활성전하, 표면 화학, 열 이력에 따라 필요한 입력이 달라진다.')
    r.p('360자리의 균일 독립 치환이라는 명시적 가정에서는1e18cm<super>-3</super>에 대해 적어도1개 도펀트가 있을 확률이 약0.7342%였다(v11). 실제 편석이나 균열확률을 측정한 값은 아니다. 이 계산 단위의 체적을 statistical correlation volume으로 자동 채택하지 않는다.','small')

    r.page('12. 열산화층 논문과 새 문헌의 범위','응력 분담, 초기 표면 반응, 기존 균열의 반응을 섞지 않았다.')
    r.p('<link href="https://doi.org/10.1016/j.engfracmech.2015.08.029" color="#087e8b">Tsuchiya 외(2016), 제공한 열산화층 논문</link>은 단조 인장 파단과 사후 파면을 다룬다. 도펀트 종·농도와 시간 분해된 최초 개시율을 주지 않는다. v11의 원문·단면·응력 분담27대조를 이번 입력 근거로 유지했다.')
    r.table(['기존 원문 감사','해석'],[
        ['bare/50/100/200nm/막 제거 강도','4.09/3.50/3.77/3.27/2.55GPa는 측정값.'],
        ['논문 FEM의 하중','측정 파단강도를 입력. 독립적인 강도 예측 검증이 아님.'],
        ['200nm와 산화막 제거 core 치수','두께를 면당/양면 합계로 읽을 때 일관되지 않는 부분이 남음.'],
        ['강도로 역산한 결함 크기','독립적으로 측정한 최초 결함 분포로 재사용하지 않음.']],[1.2,2.2])
    r.section('이번에 추가 확인한1차 자료')
    r.p('<link href="https://doi.org/10.1557/PROC-737-F8.23" color="#087e8b">Ogata 외, MRS OPL737</link>는 제목에initiation이 있지만 출판사 초록의 계산은 이미 균열이 있는 Si에서 시작한다. 초기 무균열 개시의 검증으로 채택하지 않았다.')
    r.p('<link href="https://arxiv.org/abs/0904.2091" color="#087e8b">Colombi Ciacchi 외, 습윤 산화Si 표면 연구</link>는 초기 표면 반응의 화학적 근거를 주지만 정량 반응률을 주지 않는다. <link href="https://doi.org/10.1016/j.mee.2011.04.036" color="#087e8b">Y.J.Oh 외(2012)</link>의 B 편석 연구도 개시 강도나 P/As/Sb 장벽의 대체물이 아니다.')
    r.p('PRB100,014204(2019)의 비정질silica 연구는 초록만 확인했다. 전문의 초기조건을 재현하지 못했으므로 결정Si의 개시 detector로 옮기지 않았다. 자세한 접근 범위는LITERATURE_SCOPE_V12.md에 기록했다.','small')

    r.page('13. 재현·테스트·Git 상태','계산을 다시 한 것과 저장 결과의 무결성을 확인한 것을 구분한다.')
    r.table(['항목','검증 범위'],[
        ['관련 구현 테스트',f'{validation["tests_passed"]}PASS / {validation["test_elapsed_seconds"]:.2f}s, ASE 경고2개. 최종 검증 소스 해시 포함.'],
        ['현재 소스',f'{len(source_manifest["source_files"])}개. 정적 import 추적과 실행한 기존 runner 포함.'],
        ['이번 결과',f'{len(artifact_manifest["files"])}개, {artifact_manifest["total_bytes"]/1e6:.2f}MB. 원자료·프로토콜·부분 결과 포함.'],
        ['이전 입력 자료',f'{len(source_manifest["upstream_artifact_files"])}개 해시 고정. 기존 v11 원자료와 v9 회귀 비교값.'],
        ['가중치 / 공개 자료','가중치2개와 CP2K 원문 SHA256 기록. MPA는 다운로드만 완료, 새 평가0회.'],
        ['Git 검증','작업 파일·staged blob·commit blob의 내용을 해시와 대조. 실제 결과는 최종 검증 기록을 따른다.']],[1,2.8])
    r.box('v11 source manifest는 당시 source를 기록한다. 이번에 수정한 전하 validator를 이전 해시와 같다고 주장하지 않는다. v12는 새 소스를 별도 고정하고 v11 원자료의 바이트가 유지되는지를 검증한다.')
    r.p('연구 결과 위치:results/silicon_initiation_v12<br/>재현:REPRODUCE.md 및 verify_initiation_bundle_v12.py<br/>최종 source_manifest.json/artifact_manifest.json은 실행 종료 뒤 생성했다. 계산 중 checkpoint는 완료처럼 해석하지 않는다.','small')
    r.p('외부 패키지 버전과 입력 해시를 기록했지만 완전한 운영체제·수학 라이브러리 잠금은 아니다. 원자 모델의 다른 장치·버전에서 bitwise 재현을 보장하지 않는다. PDF 옆 manifest는 실제 사용한 파일과 폰트·builder·페이지 경계를 별도로 기록한다.','small')

    r.page('14. 결론과 남은 물리적 과제','실패를 숨기지 않고 다음 계산이 필요한 이유를 좁혔다.')
    r.table(['보완한 것','이후 판단에 주는 의미'],[
        ['전체 곡률과 grip ensemble','이전 미수렴으로 알 수 없었던 특정 상태의 국소 곡률을 직접 조사할 수 있음.'],
        ['표면/내부 구조와 하중 기준','표면·내부 재배열을 구분. 기존 무하중 기준과 새100MPa 국소 접선을 분리함.'],
        ['보정 후보와 외부 전이 검사','작은 평균 오차 개선만으로 산화 계면 재료를 승인하면 안 된다는 실제 대조.'],
        ['누출 guard와 전하 경계 극한','시간 단위나 전하 평균 때문에 최초 흡수 확률이 잘못 계산되는 경로를 차단.']],[1.1,2.6])
    r.section('다음 연구의 물리 gate')
    r.p('<b>1. 재료:</b> 실제 Si/O/H 및 도펀트·전하별 intact 상태, 후보 전이상태와 끝점을 같은 전자구조 규약으로 비교한다. 공개 단일점 평균 오차만으로 장벽을 승인하지 않는다.')
    r.p('<b>2. 개시:</b> 초기 균열 없는 표면·산화층·계면에서 지속되는 새 균열 basin과 경로를 정의한다. 표면 재구성, slip, 상전이 후보와 구분하고 크기·경계·좌표 수렴을 확인한다.')
    r.p('<b>3. 동역학:</b> 유한온도 자유에너지와 대응 좌표의 실제 이동도·전이율을 얻은 뒤 같은 확률 법칙의 물리 시간을 검증한다. 그 전에는 초·Hz·피로 수명을 제시하지 않는다.')
    r.p('100MPa에서 양의 국소 곡률을 얻은 것은 endurance limit의 존재나 무한 수명을 입증한 결과가 아니다.','small')
    r.p('유한 상태의 시간 불변 generator에서 모든 초기 도달 상태가 양의율 경로로 개시 집합에 연결되면 결국 흡수된다. 진정한 비개시 집단에는 접근 불가능한 class 등 추가 조건이 필요하다. 큰 장벽이나 확률 underflow를 그 조건으로 대체하지 않는다. 자세한 가정은ENDURANCE_AND_FIRST_PASSAGE_V12.md에 기록했다.','small')
    r.box('계획한8시간 연속 실행은 배터리 최대 절전으로 완료하지 못했다. 완료한 계산·수정·검증과539행 부분 상태를 함께 제출한다. Si 웨이퍼 최초 균열의 정량 예측 모델은 아직 완성되지 않았으며 미실행 결과를 채워 넣지 않았다.')
    r.close();print(json.dumps(dict(output=args.output.name,pages=r.n,source_files=len(r.used),revision=args.revision),ensure_ascii=False))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('results','output','font','bold-font'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--as-of',required=True);p.add_argument('--label',default='최종 연구 기록')
    p.add_argument('--revision',required=True);main(p.parse_args())
