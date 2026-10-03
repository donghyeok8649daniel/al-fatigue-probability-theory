"""Korean report built from verified raw tables, with no new scientific fits."""
import argparse,json,os
from pathlib import Path
from xml.sax.saxutils import escape
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak,Image


def main(args):
    root=args.package;read=lambda name:json.loads((root/name).read_text(encoding='utf-8'))
    if args.font is None or args.bold_font is None:
        system_root=os.environ.get('SystemRoot')
        if not system_root:raise ValueError('specify --font and --bold-font with Korean TTF files')
        args.font=args.font or Path(system_root)/'Fonts/malgun.ttf'
        args.bold_font=args.bold_font or Path(system_root)/'Fonts/malgunbd.ttf'
    pdfmetrics.registerFont(TTFont('Korean',str(args.font)))
    pdfmetrics.registerFont(TTFont('KoreanBold',str(args.bold_font)))
    pdfmetrics.registerFontFamily('Korean',normal='Korean',bold='KoreanBold',italic='Korean',boldItalic='KoreanBold')
    ink=colors.HexColor('#172c3b');accent=colors.HexColor('#146580');muted=colors.HexColor('#526777')
    styles={
        'body':ParagraphStyle('body',fontName='Korean',fontSize=10.4,leading=16.7,textColor=ink,spaceAfter=8,wordWrap='CJK'),
        'small':ParagraphStyle('small',fontName='Korean',fontSize=8.3,leading=12.7,textColor=muted,spaceAfter=6,wordWrap='CJK'),
        'title':ParagraphStyle('title',fontName='KoreanBold',fontSize=25,leading=34,textColor=ink,spaceAfter=13,wordWrap='CJK'),
        'heading':ParagraphStyle('heading',fontName='KoreanBold',fontSize=17,leading=24,textColor=accent,spaceAfter=13,wordWrap='CJK'),
        'sub':ParagraphStyle('sub',fontName='KoreanBold',fontSize=11.4,leading=18,textColor=ink,spaceAfter=7,spaceBefore=8,wordWrap='CJK'),
        'cell':ParagraphStyle('cell',fontName='Korean',fontSize=8.2,leading=12,textColor=ink,wordWrap='CJK'),
    }
    story=[]
    def p(text,style='body'):story.append(Paragraph(text,styles[style]))
    def h(text):p(text,'heading')
    def sub(text):p(text,'sub')
    def page():story.append(PageBreak())
    def image(name,width=485):
        im=Image(str(root/'figures'/name));im.drawHeight=im.imageHeight*width/im.imageWidth;im.drawWidth=width
        story.append(im);story.append(Spacer(1,7))
    def table(rows,widths,padding=7):
        data=[[Paragraph(escape(str(cell)),styles['cell']) for cell in row] for row in rows]
        t=Table(data,colWidths=widths,hAlign='LEFT',repeatRows=1)
        t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e4eef2')),
            ('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),7),('RIGHTPADDING',(0,0),(-1,-1),7),
            ('TOPPADDING',(0,0),(-1,-1),padding),('BOTTOMPADDING',(0,0),(-1,-1),padding),
            ('LINEBELOW',(0,0),(-1,0),.7,accent),('LINEBELOW',(0,1),(-1,-1),.25,colors.HexColor('#cdd8df'))]))
        story.append(t);story.append(Spacer(1,10))
    selection=read('final_selection.json')
    metrics=read('reaggregated_group_metrics.json');selected=read('raw_runs/'+selection['selected_run']+'/selected_fit.json')
    phonon=read('raw_runs/final_selected_phonon_recovery_results/summary.json');verification=read('independent_verification.json')
    def m(case,kind,key):return next(r[key] for r in metrics if r['case']==case and r['config_type']==kind and r['declared_xc']=='PW91')
    a=m('original_SW','surface_111','force_RMSE_eV_A');b=m('joint_static_candidate','surface_111','force_RMSE_eV_A')
    ea=m('original_SW','surface_111','relative_energy_RMSE_eV_atom');eb=m('joint_static_candidate','surface_111','relative_energy_RMSE_eV_atom')
    p('실리콘 웨이퍼 연구 현황','title')
    p('원자 에너지의 베셀 표현과 물성 보정','heading')
    p(args.as_of+' 기준 · 계산 결과와 독립 재검증을 함께 정리','small')
    sub('현재 어디까지 왔나')
    p('순수 실리콘의 Stillinger-Weber(SW) 에너지를 푸리에-베셀 합으로 계산하는 구현을 끝냈다. 같은 원자 배열에서 직접 합과 베셀 합의 에너지, 힘, 전체 곡률을 대조했다. 이번 작업에서는 이 표현을 유지하면서 실리콘의 벌크 탄성, 표면 힘, 상대 에너지 오차를 줄이는 보정을 진행했다.')
    p('현재 후보는 평형 격자와 세 탄성값을 맞추고, 벌크 상대 에너지 오차를 원래 SW보다 키우지 않는 조건을 동시에 만족한다. 앞선 후보들에서 충돌하던 정적 조건을 함께 맞춘 결과다. 학습에 넣지 않은 Si(111) 표면에서 일부 개선을 얻었지만, 다른 구조와 진동수에는 오차가 남아 재료 모델로 채택하지 않았다.')
    table([['검사 항목','원래 SW','현재 비교 후보'],
        ['Si(111) 힘 RMSE',f'{a:.4f} eV/Å',f'{b:.4f} eV/Å · {(1-b/a)*100:.1f}% 감소'],
        ['Si(111) 상대 에너지 RMSE',f'{ea:.4f} eV/원자',f'{eb:.4f} eV/원자 · {(1-eb/ea)*100:.1f}% 감소'],
        ['벌크 상대 에너지 RMSE','0.02881 eV/원자',f"{m('joint_static_candidate','dia','relative_energy_RMSE_eV_atom'):.5f} eV/원자"],
        ['벌크 힘 RMSE','0.17690 eV/Å',f"{m('joint_static_candidate','dia','force_RMSE_eV_A'):.5f} eV/Å · 악화"]], [197,135,155])
    p('실제 첫 균열 개시 확률, 피로수명, 물리적 초·Hz는 아직 확보하지 못했다. 이번 보정은 기존 DFT 자료와 대조한 정적 계산이다. 새 DFT나 MD 궤적은 만들지 않았다.','small')
    sub('이 보고서에서 구분하는 것')
    p('계산식의 정확성, 재료 물성의 정확성, 균열 개시 동역학의 정확성은 각각 확인해야 한다. 베셀 구현과 회귀시험 통과는 첫 번째에 해당한다. 이번 작업은 두 번째를 개선하고 그 한계를 확인한 단계다.')
    page()
    h('1. 지배방정식과 실리콘 에너지의 역할')
    p('알루미늄 연구에서 쓰던 확률밀도 진화의 구조를 유지한다. 상태 q의 에너지와 열운동이 확률 흐름을 만들고, 개시 경계로 빠져나가는 질량을 누적한다. 원자 에너지나 소재만 바꿨다고 이 전체 모델이 검증되는 것은 아니다.')
    table([['공통 확률식','의미'],
        ['∂t P = -Σi ∂qi Ji','상태공간에서 확률밀도가 이동한다.'],
        ['Ji = -Σj Mij [P ∂qj G + kBT ∂qj P]','에너지 기울기의 구동과 열확산을 함께 사용한다.'],
        ['G = F(q,T) - 외력의 일','F는 선택한 집단좌표의 자유에너지다.'],
        ['Pinit = 누적 개시 경계 흡수 질량','수치 질량 오차와 물리적 흡수를 분리한다.']], [258,229])
    p('현재 Si 정적 계산의 좌표는 q=(δa, u1, u2)다. δa는 한 (111) 계면의 벌어짐, u1과 u2는 두 횡방향 이동이다. 두 반결정을 강체로 두고, 원자마다 모든 이웃의 상호작용을 합한 뒤 같은 벌크 상태를 빼서 계면 초과 에너지 W(q)를 구한다.')
    p('법선·전단 하중을 이 좌표에 걸 때 정적 일은 A0(Tn δa + τ1 u1 + τ2 u2)다. A0는 명시된 원자 계면 셀 면적이다. 시편 전체 응력, 이 좌표의 하중, 집단좌표의 에너지 정규화를 먼저 맞춰야 한다. 통계적 상관 면적 Ac는 이 에너지식에 들어가지 않는다.')
    sub('kBT를 고체에 쓸 수 있나')
    p('고체도 열평형에서 열적 요동을 갖는다. kBT를 쓸 수 있는지는 고체의 강성 자체보다 선택한 좌표와 소거한 자유도의 통계에 달려 있다. 해당 좌표의 자유에너지와 이동도, 필요한 시간 분리를 확인해야 한다. 원자들의 강한 상관과 빠른 진동을 소거한 뒤에도 기억효과가 남으면, 단순한 일정 이동도 확률식으로 충분한지 추가 검증이 필요하다.')
    p('이번에 계산한 W는 0 K 정적 에너지다. 유한 온도에서 다른 원자 자유도를 적분한 F(q,T), 이동도 M, 물리 시간 환산은 아직 구하지 않았다. 조화 진동수의 THz를 확률식의 이동도나 피로시험 Hz로 대입하지 않는다.','small')
    p('여러 원자의 집단좌표 에너지를 원자당으로 나누면, 같은 확률법칙을 유지하도록 열항과 이동도도 함께 환산해야 한다. 에너지 척도만 바꾸고 kBT를 그대로 유지한 결과를 같은 동역학으로 보지 않는다.','small')
    page()
    h('2. 베셀 표현에서 실제로 구현한 것')
    p('Si는 다이아몬드 결정 구조를 가지며, SW는 원자 두 개의 거리항과 중심 원자-이웃 두 개의 결합각 항을 갖는다. Al의 LJ 식을 그대로 옮기면 이 결합각 정보를 잃는다. 현재 구현은 SW의 결합각 항을 유지한 채 주기적인 면 합을 변환한다.')
    table([['에너지 구성','현재 계산'],
        ['거리항','원래 SW의 유한 지지 함수와 거리 미분을 유지한다.'],
        ['각도항','g(rij)g(rik)(cos θjik - c0)²를 중심 원자마다 합한다.'],
        ['각도 형상 확장','x=cos θ-c0, f(x)=x²(1-βx)². β=0이면 원래 SW다.'],
        ['베셀 적분','각도 차수에 맞춘 J0부터 J4까지의 Hankel 적분과 역격자 합.']], [139,348])
    p('모든 이웃 면에서 거리 가중치와 방향 텐서를 먼저 합하고, 그 다음 텐서 제곱을 계산한다. 같은 이웃을 두 번 세는 자기항을 빼고, 각 중심 원자의 벌크 항을 차감한다. 다이아몬드 구조의 3차 방향 텐서는 0이 아니므로 이를 남겼다.')
    p('확장 각도항은 cos θ의 4차 다항식이다. 대칭 텐서의 고정 성분 35개로 정확히 표현할 수 있다. 이 35개는 독립 적합 계수가 아니라, β 한 개를 계산하는 통계량이다. 같은 면의 비자기 이웃은 직접 합하고, 면 간 항은 베셀 합으로 계산한다.')
    sub('수치 검사')
    p('기존 SW 베셀 구현은 12개 상태에서 직접 합과 대조했다. 이번 후보와 추가 이웃·선호각 변경에는 별도의 직접 합, 보존적 힘, 변형률 미분, 전체 Hessian, 베셀 대조를 수행했다. 곡률 행렬은 법선과 두 횡방향의 교차항까지 포함한다.')
    p('이번 각도 후보 대조의 허용 범위는 에너지 2×10^-8 eV, 기울기 2×10^-7 eV/Å, Hessian 3×10^-6 eV/Å²다. 검사한 상태에서 통과했다. 역격자 셸과 구적의 유한 분해능을 사용하므로, 이 결과를 모든 상태의 보편 오차 보증으로 확대하지 않는다.','small')
    resolution=read('raw_runs/final_selected_convergence_results/summary.json')
    p(f"최종 후보에서 역격자 셸 48/96/128과 구적 256/512/768을 따로 바꿔 20개 베셀 계산을 추가했다. 법선 임계 상태와 횡이동 상태를 직접 합과 대조했다. 가장 세밀한 해상도의 최대 곡률 차이는 {resolution['final_refinement_maximum_errors']['hessian_error_eV_A2']:.2e} eV/Å²다. 해상도별 오차가 완전히 단조롭지는 않으며, 전 상태의 꼬리 오차 상한을 인증한 것은 아니다.",'small')
    page()
    h('3. 보정에 쓴 자료와 선택 규칙')
    p('공개된 Si GAP 구축용 원자료 2,475개 구조, 171,815개 원자를 사용했다. GAP를 새로 학습한 작업은 아니다. 구조에 기록된 교환-상관 함수 PW91, PBE, 미상 자료를 분리해 다뤘다.')
    table([['역할','자료'],
        ['학습 511개','PW91 벌크 489개, (001) 표면 10개, (110) 표면 12개.'],
        ['손실에서 제외 1,964개','(111) 표면 전체, 다른 결정상·결함·분리·기존 균열을 포함한 나머지 구조.'],
        ['알려진 벌크 보정 기준','공개 테스트 자료의 격자, C11/C12/C44, 완전 결정 E(V).'],
        ['선택 이후 외부 검사','별도 공개 세트 28개 구조, 1,020개 원자. 선택과 손실에 사용하지 않았다.']], [140,347])
    p('힘·상대 에너지·정적 응력을 공동 손실에 넣고, 학습 구조군마다 같은 총 가중치를 부여했다. 상대 에너지는 같은 PW91 벌크 구조 하나를 고정 기준으로 썼다. 이 구조는 열적으로 변형된 프레임이며, 완전 이완한 기저상태나 응집 에너지 기준이 아니다.')
    p('세 탄성 기준은 공개된 유한 변형 응력 적합값이다. 정확한 DFT 0변형 접선이 독립 인증된 값과는 구별한다.','small')
    p('손실의 척도는 힘 0.2 eV/Å, 상대 에너지 0.05 eV/원자, 응력 0.05 eV/Å³, 완전 결정 E(V) 0.02 eV/원자다. 선택을 위한 선언된 가중치이며, DFT의 불확실성이나 통계적 신뢰구간을 뜻하지 않는다.')
    sub('벌크를 악화시키지 않는 조건')
    p('벌크 상대 에너지 오차는 원래 SW의 RMSE 0.0288054 eV/원자를 넘지 않게 제한했다. 거리·각도 진폭은 양수로 유지하고, 내부 이완과 탄성의 안정성도 검사했다. 조건을 통과한 후보끼리 학습 손실로 선택했다. 제외 구조에서 좋아 보이는 후보를 다시 골라 최종 결과로 만들지 않았다.')
    p('제외 자료와 외부 28개는 앞선 후보에서 이미 검사한 자료다. 최종 후보의 손실·선택에는 쓰지 않았지만, 완전히 새 독립 실험은 아니다. 외부 자료의 정렬한 주기 좌표와 셀 지문에서 학습 자료와 중복이 없었으나, 모든 회전·평행이동·대체 셀의 동등성까지 증명한 것은 아니다.','small')
    page()
    h('4. 어떤 보정이 개선됐고 무엇이 실패했나')
    labels=[('original_SW','원래 SW'),('amplitude_fit','진폭만 보정'),('bulk_energy_guard','벌크 에너지 제한'),
            ('angular_grid','각도 격자 후보'),('angular_shape','앞선 형상 · C44 부족'),
            ('joint_static_candidate','현재 공동 탄성 후보'),
            ('extended_support_diagnostic','추가 이웃 · 제한 실패'),('angle_preference_diagnostic','선호각 변경 · 제한 실패')]
    rows=[['후보','벌크 F','벌크 ΔE','(111) F','(111) ΔE']]
    for case,label in labels:rows.append([label,f"{m(case,'dia','force_RMSE_eV_A'):.3f}",f"{m(case,'dia','relative_energy_RMSE_eV_atom'):.4f}",
        f"{m(case,'surface_111','force_RMSE_eV_A'):.3f}",f"{m(case,'surface_111','relative_energy_RMSE_eV_atom'):.4f}"])
    table(rows,[203,69,75,69,71],padding=5)
    p('F는 힘 RMSE(eV/Å), ΔE는 고정 기준에 대한 상대 에너지 RMSE(eV/원자)다. 표면 힘이 크게 줄어든 첫 진폭 보정은 벌크를 악화시켰다. 벌크 에너지를 지킨 보정도 탄성을 악화시켰다.','small')
    image('material_errors.png')
    p(f"현재 후보는 σ={selected['parameters']['sigma']} Å, γ={selected['parameters']['gamma']}, c0={selected['parameters']['costheta0']:.6f}, β={selected['angle_beta']:.6f}다. 격자 5.461022 Å와 C11/C12/C44=153.28991/56.25009/72.17693 GPa를 맞췄다. 바꾼 선호각의 비영 벌크 항까지 다시 계산했다. 앞선 β=0.375 후보는 C44가 9.6% 낮지만 표면 오차는 더 작다.")
    p(f"추가 이웃 16개, 선호각 단독 9개, 초기 공동 19개는 벌크 에너지 제한에 실패했다. 양의 각도 경계와 주변 형상 {selection['examined_profiles']}개 중 {selection['eligible_profiles']}개가 선언한 조건을 통과했고 학습 손실로 선택했다. 전역 최적성이나 모든 환경의 정확성을 뜻하지 않는다.",'small')
    p('선택 후 같은 세 힘 기저를 고정한 진단에서는, 탄성·에너지 조건을 풀고 Si(111) 힘만 맞추면 RMSE가 1.021에서 0.489 eV/Å로 내려갔다. 이 계수는 재료 조건을 유지한 적합이 아니다. 표면 오차의 일부가 공동 물성 제약과 충돌함을 보여주며, 다른 형상이나 전체 에너지 계열의 불가능성 증명은 아니다.','small')
    p('저장된 벌크 행렬로 따로 대조하면, 표면 힘만 맞춘 계수는 기준 격자에서 법선 응력 -1.501 GPa를 남겨 무하중 평형을 벗어난다. 현재 형상에서는 평형과 두 법선 곡률의 세 독립 조건이 세 진폭을 정한다. 진폭만 더 조정하면 이 조건을 함께 유지할 수 없다.','small')
    page()
    h('5. 외부 구조와 진동수에서 드러난 한계')
    image('external_errors.png')
    p('선택 이후의 외부 세트에서는 두 침입형 원자 결함의 힘 오차가 줄었다. 입계, 비정질, 적층결함은 악화됐다. 구조 종류에 따라 개선 방향이 다르므로, 현재 후보가 Si의 여러 환경을 함께 정확히 설명한다고 보기 어렵다.')
    image('phonon_bands.png')
    p(f"최종 선택 후보를 새로 조립한 조화 진동수의 DFT 대비 RMSE는 {phonon['band_RMSE_THz']:.3f} THz, 최대 차이는 {phonon['band_maximum_absolute_error_THz']:.3f} THz다. Γ 광학 모드는 {phonon['gamma_optical_THz'][0]:.3f} THz이며 DFT {phonon['reference_gamma_optical_THz'][0]:.3f} THz와 차이가 남는다.")
    p('문헌은 4×4×4 초셀의 0.03 Å 유한 변위, 현재 계산은 정확한 국소 Hessian을 사용했다. 공개 SW 대조와의 밴드 차이는 최대 0.00477 THz여서, 현재 DFT 차이를 구현의 단순 단위 오류로 설명할 수는 없다. 9³ 파수 표본 안정성을 검사했으나 전체 Brillouin 영역의 안정성 증명은 아니다.','small')
    page()
    h('6. 계면 곡선과 첫 균열 개시는 구별했다')
    image('rigid_interface.png')
    branch=[r for r in read('raw_runs/interface_recovery_results/branch_roots.json') if r['case']=='original_SW']+read('raw_runs/final_selected_interface_results/branch_roots.json');rows=[['강체 계면','최초 국소 안정성 상실 하중','강체 분리 일']]
    for row in branch:
        if row['mesh_points']!=161:continue
        rows.append([('원래 SW' if row['case']=='original_SW' else '각도 후보')+' / '+row['cut_kind'],
            f"{row['normal_traction_GPa']:.3f} GPa",f"{row['rigid_separation_work_J_m2']:.3f} J/m²"])
    table(rows,[155,192,140])
    p('닫힌 완전 결정 상태에서 법선 하중을 늘린 강체 대칭 경로를 검사했다. 법선과 두 횡방향의 3×3 Hessian 최소 고유값이 처음 0이 되는 구간을 추적했다. 이 경로에서는 처음 상실하는 모드가 법선 방향이었다.')
    p('이 값은 정해진 원자 평면을 두 강체 반결정으로 분리할 때의 국소 기계적 값이다. 표면 원자 이완, 국소 균열핵, 열적 전이, 시편 크기와 결함 분포가 없다. 실제 첫 균열 개시 강도나 수명으로 사용하지 않는다. MPa 영역의 시편 거동을 설명하려면 이런 상태와 시편 하중을 연결하는 원자료가 필요하다.')
    p('첫 표시에서 응력 변환에 표면 에너지 단위를 잘못 적용한 오류를 발견했다. 원래 에너지·기울기·Hessian을 보존하고, SI 단위에서 다시 계산했다. 1 eV/Å³=160.2176634 GPa, 1 eV/Å²=16.02176634 J/m²를 분리했다. 이 페이지와 최종 교정 자료의 하중은 올바른 값이다.','small')
    page()
    h('7. 수치 신뢰성과 미해결 참조 차이')
    git_tests=read('tests_receipt.json')['git_index_replay_run']
    table([['검증','실제 결과'],
        ['실제 Git 바이트의 신규 전체 시험',f"{git_tests['passed']}개 통과 · {git_tests['elapsed_seconds']:.2f}초."],
        ['새 보정·각도·진동수 검사','28개 통과 · 163.56초.'],
        ['선호각·각도 형상의 비영 벌크 항 검사','5개 통과 · 122.88초.'],
        ['기존 환경·기하·격자 회귀','24개 통과 · 6.87초. 별도 실행.'],
        ['보존된 원자 예측 재집계',f"2,475개 구조 · {verification['raw_prediction_cases_verified']}개 후보 · {verification['group_metrics_verified']}개 구조군 일치."],
        ['과거 계산 소스 검사',f"{verification['exact_source_bindings_verified']}개 해시 연결 확인."],
        ['축약 적합 손실 재생',f"{verification['compact_Gram_profiles_replayed']}개 프로파일의 Gram 통계 재생."],
        ['현재 계면 곡률 독립 차분',f"최대 차이 {verification['new_candidate_force_Hessian_maximum_difference_eV_A2']:.2e} eV/Å²."]], [219,268])
    p('새 힘 기저로 원래 SW를 복원한 결과는 기존 독립 LAMMPS 원자료와 최대 힘 차이 3.55×10^-13 eV/Å, 에너지 차이 1.52×10^-11 eV였다. 이는 같은 SW 계산을 재현했다는 검사다. DFT 물성 오차가 해결됐다는 뜻은 아니다.')
    fits=read('raw_runs/final_selected_finite_strain_results/fits.json');fine=next(r for r in fits if r['pattern']=='cubic_exx_plus_yz' and r['increment']==.001)
    p(f"현재 후보의 내부 원자힘을 따로 이완한 응력-변형률 검사에서도 세 탄성값이 일치했다. 변형 간격을 1%에서 0.1%로 줄이면 C11/C12/C44는 {fine['C11_GPa']:.4f}/{fine['C12_GPa']:.4f}/{fine['C44_GPa']:.4f} GPa로 접선값에 접근한다.",'small')
    sub('계산 중 발견해 고친 문제')
    p('압축에서 두 번째 이웃이 들어오면 “벌크 곡률 보존” 가정이 깨졌다. 완전 결정 E(V) 재사용을 감사해 새로 864점을 계산했고, 이전 선택은 그대로였다. 각도 연속 적합의 E(V) 정규화 오류도 별도 교정했다. 마지막 JSON 저장 실패가 난 형상 계산은 원자료·선택 잠금·소스 해시로 복구했고, 잃어버린 실행 시간은 미상으로 남겼다.')
    sub('문헌 SW의 전단계수는 아직 완전히 재현하지 못했다')
    p('문헌의 결합 변형 패턴을 따라 내부 원자힘을 충분히 이완하면 C44는 56.323 GPa(1% 간격), 56.427 GPa(0.1% 간격)였다. 0변형 접선은 56.429 GPa다. 문헌의 53.367 GPa는 이 변형 폭만으로 재현되지 않았다. 과거 quippy와 최적화 환경 전체를 복원한 결과는 없으며, 이 차이로 DFT 탄성값을 임의 교정하지 않았다.')
    page()
    h('8. 첫 균열 개시 모델로 이어가기 위해 남은 것')
    p('목표는 미리 균열을 넣은 시편의 성장 속도가 아니라, 초기 시편에서 처음 지속되는 균열핵이 생기는 조건과 확률을 구하는 것이다. 한 결합의 순간적 늘어남, 표면 재구성, 전체 평면의 분리 장벽을 같은 개시 사건으로 세지 않는다.')
    table([['필요한 근거','현재 상태'],
        ['초기 상태 A와 첫 개시 상태 B','공간적 연결, 지속성, 복귀 가능성을 구별하는 원자 사건 정의와 독립 경로 검증이 더 필요하다.'],
        ['관련 환경의 재료 에너지','벌크·계면·표면·결함을 함께 만족하는 후보를 확보하지 못했다.'],
        ['자유에너지 F(q,T)','정적 W와 구별한 집단좌표 자유에너지 계산이 필요하다.'],
        ['이동도 M과 물리 시간','좌표에 대응하는 응답·상관·기억효과 검증이 필요하다. 현재 미확보.'],
        ['도핑·산화층·환경','이번 적합은 순수 Si다. 도펀트 물질·농도·전하, 산화층 및 수분·잔류응력의 새 물성 보정은 없다.'],
        ['시편 확률과 실험 검증','공간 하중, 상관 면적, 개시 자료의 측정 조건과 표본수가 필요하다.']], [163,324])
    sub('이번 결과로 다음에 판단할 수 있는 것')
    p('벌크 에너지를 지키지 않은 낮은 표면 손실과, 탄성값만 맞춘 후보를 빠르게 걸러낼 수 있게 됐다. 베셀 합과 직접 원자힘을 같은 에너지에서 대조할 기반도 생겼다. 다음 재료 후보는 이 조건들과 파수별 진동수까지 함께 검사해야 한다.')
    p('재료 후보가 이 검사를 통과한 뒤에, 초기 무균열 구조에서 국소 개시 상태와 전이 경로를 확보해야 한다. 그 경로를 바탕으로 자유에너지와 동역학을 계산하고, 마지막에 시편 조건과 비교한다. 현재 비교 후보의 정적 결과를 먼저 생산 확률식에 넣어 수명을 출력할 근거는 없다.','small')
    page()
    h('9. 원자료, 출처와 재현 방법')
    p('산출물에는 선택 후보의 원자별 예측, 구조군별 오차, 각 단계의 작은 원자료, 당시 계산 소스와 SHA-256, 적합 손실을 재생할 Gram 통계, 독립 검증기를 보존했다. 큰 중복 설계 배열은 기존 작업 캐시에 남겼고, 전체 목록·크기·해시를 함께 기록했다.')
    table([['파일','확인할 내용'],
        ['selected_predictions.npz',f"{verification['raw_prediction_cases_verified']}개 후보의 2,475개 구조 에너지와 원자별 힘, DFT 참조."],
        ['reaggregated_group_metrics.json','PW91/PBE/미상 및 학습 여부를 나눈 오차.'],
        ['raw_runs/','실제 계산·실패·교정 결과와 원래 소스 스냅샷.'],
        ['design_gram_replay.json','형상별 손실·벌크 제한을 원자료 통계로 재생.'],
        ['verify_material_package.py','해시·원자 예측·손실·계면 SI 단위를 독립 검사.'],
        ['README.md','의존성, 실행 인자와 새 적합을 다시 수행하는 순서.']], [206,281])
    sub('자료 출처')
    for text in [
        'Bartók 외, Physical Review X 8, 041048 (2018). Si 원자료와 검증 체계. DOI: <link href="https://doi.org/10.1103/PhysRevX.8.041048" color="#146580">10.1103/PhysRevX.8.041048</link>.',
        'Cambridge 공개 원자 자료. DOI: <link href="https://doi.org/10.17863/CAM.65004" color="#146580">10.17863/CAM.65004</link>. 기록된 구조·교환-상관 함수 정보를 사용했다.',
        'libAtoms silicon-testing-framework의 고정 공개 태그. DOI: <link href="https://doi.org/10.5281/zenodo.1250555" color="#146580">10.5281/zenodo.1250555</link>. 고정 Git 트리 fc252cb7d41df7e2bc672d614f3d76a40c9f2ecb.',
        'SW 원래 순수 Si 파라미터는 기존 저장소의 source_Si.sw. 이번의 원래 SW ε=2.1683 eV와 공개 프레임워크 ε=2.1675 eV를 구별했다.',
        '과거 matscipy 탄성 코드의 참조 커밋은 1202c4e9c85799ee15028faf7e7e4b4d0436fe72. 실제 문헌 실행 환경이 이 커밋으로 고정됐다는 증거는 없다.',
    ]:p(text,'small')
    p('이 보고서의 개선율은 위 원자료에서 계산했다. 물성 채택, 첫 개시 검증, 물리 시간 보정은 모두 미승인 상태로 보존했다. Al·UI·생산 에너지와 시간 설정은 이번 연구에서 변경하지 않았다.','small')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    def footer(canvas,doc):
        canvas.saveState();w,h=A4;canvas.setStrokeColor(accent);canvas.setLineWidth(.8);canvas.line(52,h-40,w-52,h-40)
        canvas.setFont('Korean',8);canvas.setFillColor(muted);canvas.drawString(52,h-31,'실리콘 웨이퍼 · 원자 에너지와 물성 보정')
        canvas.drawString(52,28,args.as_of);canvas.drawRightString(w-52,28,str(doc.page));canvas.restoreState()
    doc=SimpleDocTemplate(str(args.output),pagesize=A4,rightMargin=52,leftMargin=52,topMargin=57,bottomMargin=48,
        title='실리콘 웨이퍼 연구 현황',author='실리콘 웨이퍼 연구',pageCompression=1)
    doc.build(story,onFirstPage=footer,onLaterPages=footer)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--package',type=Path,default=Path(__file__).resolve().parent);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--as-of',default='2026년 10월 3일')
    parser.add_argument('--font',type=Path)
    parser.add_argument('--bold-font',type=Path)
    main(parser.parse_args())
