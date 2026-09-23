"""Build the Korean v11 research report from saved evidence; no new physics run.

ReportLab and explicitly supplied Korean fonts are required. Page layout is
bounded and image/text pages must still be rendered and visually reviewed.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path
from xml.sax.saxutils import escape
from reportlab.pdfgen import canvas
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import Paragraph,Table,TableStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.utils import ImageReader

W,H=595.276,841.89
L,CW=46,503.276
INK=colors.HexColor('#193448');TEAL=colors.HexColor('#087e8b')
GRAY=colors.HexColor('#506674');LINE=colors.HexColor('#d9e3e8')


class Report:
    def __init__(self,args):
        self.args=args;self.used={};self.pages=[];self.n=0;self.y=0
        pdfmetrics.registerFont(TTFont('Korean',str(args.font)))
        pdfmetrics.registerFont(TTFont('KoreanBold',str(args.bold_font)))
        pdfmetrics.registerFontFamily('Korean',normal='Korean',bold='KoreanBold',italic='Korean',boldItalic='KoreanBold')
        self.styles={k:ParagraphStyle(k,fontName='KoreanBold' if k in ('title','section','head') else 'Korean',
            fontSize=size,leading=leading,textColor=colors.white if k=='head' else INK,wordWrap='CJK')
            for k,size,leading in [('title',20,28),('section',12,18),('body',10,16),('small',8.5,13),('cell',8.5,12.5),('head',8.5,12.5)]}
        args.output.parent.mkdir(parents=True,exist_ok=True)
        self.c=canvas.Canvas(str(args.output),pagesize=(W,H),pageCompression=1)
        self.c.setTitle('Si 웨이퍼 균열 개시 연구 v11 - 계산·검증·한계')
        self.c.setAuthor('Si wafer research project')
    def use(self,path):
        path=Path(path);relative=path.relative_to(self.args.results).as_posix()
        self.used[relative]=dict(bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        return path
    def data(self,relative):return json.loads(self.use(self.args.results/relative).read_text(encoding='utf-8'))
    def page(self,title,subtitle):
        if self.n:self.finish();self.c.showPage()
        self.n+=1;self.title=title
        self.c.bookmarkPage(f'p{self.n}');self.c.addOutlineEntry(title,f'p{self.n}',0,False)
        self.c.setFillColor(TEAL);self.c.rect(L,H-41,33,3,fill=1,stroke=0)
        self.c.setFont('KoreanBold',8);self.c.drawString(L+43,H-43,'SILICON / CRACK INITIATION / v11')
        self.y=H-68;self.p(title,'title',8);self.p(subtitle,'small',18)
    def finish(self):
        if self.y<65:raise ValueError(f'page {self.n} overflow: y={self.y}, {self.title}')
        self.pages.append(dict(page=self.n,title=self.title,content_bottom=self.y))
        self.c.setStrokeColor(LINE);self.c.line(L,51,W-L,51)
        self.c.setFillColor(GRAY);self.c.setFont('Korean',7)
        self.c.drawString(L,38,f'{self.args.as_of} KST | {self.args.label}')
        self.c.drawRightString(W-L,38,f'{self.n:02d} / 16')
    def p(self,text,style='body',after=9):
        p=Paragraph(text,self.styles[style]);_,height=p.wrap(CW,1000)
        p.drawOn(self.c,L,self.y-height);self.y-=height+after
    def section(self,text):self.y-=3;self.p(text,'section',7)
    def box(self,text):
        p=Paragraph(text,self.styles['body']);_,height=p.wrap(CW-24,1000)
        self.c.setFillColor(colors.HexColor('#edf5f7'));self.c.roundRect(L,self.y-height-20,CW,height+20,5,fill=1,stroke=0)
        p.drawOn(self.c,L+12,self.y-height-10);self.y-=height+32
    def table(self,head,rows,widths):
        cells=[[Paragraph(str(x),self.styles['head']) for x in head]]
        cells += [[Paragraph(str(x),self.styles['cell']) for x in row] for row in rows]
        table=Table(cells,colWidths=[CW*w/sum(widths) for w in widths])
        table.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),TEAL),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f0f5f7')]),
            ('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),7),('RIGHTPADDING',(0,0),(-1,-1),7),
            ('TOPPADDING',(0,0),(-1,-1),7),('BOTTOMPADDING',(0,0),(-1,-1),7),('LINEBELOW',(0,0),(-1,0),.5,TEAL)]))
        _,height=table.wrap(CW,1000);table.drawOn(self.c,L,self.y-height);self.y-=height+13
    def picture(self,relative,max_height=245):
        path=self.use(self.args.results/relative);img=ImageReader(str(path));w,h=img.getSize()
        scale=min(CW/w,max_height/h);dw,dh=w*scale,h*scale
        self.c.drawImage(img,L+(CW-dw)/2,self.y-dh,width=dw,height=dh,mask='auto');self.y-=dh+10
    def close(self):
        self.finish()
        if self.n!=16:raise ValueError('unexpected page count')
        self.c.save()
        manifest=dict(as_of_KST=self.args.as_of,label=self.args.label,pages=self.pages,
            output_filename=self.args.output.name,pdf_sha256=hashlib.sha256(self.args.output.read_bytes()).hexdigest(),source_files=self.used,
            source_revision=self.args.revision,visual_QA_required=True,
            builder_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            font_sha256=hashlib.sha256(self.args.font.read_bytes()).hexdigest(),
            bold_font_sha256=hashlib.sha256(self.args.bold_font.read_bytes()).hexdigest())
        self.args.output.with_suffix('.manifest.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')


def main(args):
    r=Report(args);root=args.results
    source=r.data('thermal_oxide_source.json');surface=r.data('oxidized_surfaces_max320/summary.json')
    env=r.data('surface_environment_summary_v2/summary.json');energy=r.data('surface_energy_differences/summary.json')
    force=r.data(args.force_replay+'/summary.json');math=r.data('first_passage_math/summary.json')
    occupancy=r.data('doping_occupancy/summary.json');validation=r.data('validation.json')
    displacement=[]
    for path in sorted((root/'intact_prism_360_linesearch_v2').glob('state_*/result.json')):
        displacement.append(json.loads(r.use(path).read_text(encoding='utf-8')))
    converged=[row for row in displacement if row['converged']]
    continuation=r.data('prism_10pct_continuation/summary.json')
    continuation_audit=r.data('prism_10pct_replay/summary.json')
    converged_strains={round(row['strain'],12) for row in converged}
    if continuation['converged']:converged_strains.add(round(continuation['strain'],12))
    qecurve=r.data('oxide_mace_mtpu_1159_v2/summary.json')
    surface_group=next(g for g in surface['groups'] if g['group']=='Si_O_H_vacuum1')
    r.page('Si 웨이퍼 균열 개시 연구','v11 | 열산화층 논문 검토에서 원자 기준 계산과 공통 확률 구조까지')
    r.box('<b>목표는 미리 균열을 넣지 않은 Si 초기 시편의 첫 균열 개시다.</b><br/>현재 확보한 것은 원문·재료 에너지·수치 구현·무축하중 기준에 대한 검증이다. 실제 Si의 개시 장벽, 개시확률, 물리 수명은 아직 검증되지 않았다.')
    r.table(['이번에 완료한 일','직접 확인한 결과'],[
        ['열산화층 원문 감사','5개 실험그룹, 3개 기하규약 × 3개 잔류응력 해석. 측정 강도를 FEM 하중으로 재사용했음을 확인.'],
        ['외부 DFT와 실제 MACE 대조','QE 1,159상태 전체, CP2K 1,013상태/179,857원자. 독립 parser로 raw 배열과 지표 재검증.'],
        ['균열 seed 없는 360원자 시편',f'후속 이완을 포함해 수렴 strain {len(converged_strains)}개. 별도 축력0/100MPa 두 상태 완료.'],
        ['원인별 첫 통과 확률 구현',f'관련 {validation["tests_passed"]}개 테스트 통과. 합성 generator에서만 검증했으며 Si rate는 미보정.']],[1,2.45])
    r.section('가장 중요한 결론')
    r.p('산화 표면의 산소 이웃1~3개 Si에서 큰 힘 오차가 유지된다. 같은 조성의 에너지 순서도 일부 틀린다. 따라서 현재 중성 MACE의 결과를 바로 개시 장벽이나 실제 웨이퍼 수명으로 승격할 근거가 부족하다.')
    r.p('Al/UI/생산 에너지와 physical-Hz gate는 보존했다. 새 DFT·MD·재학습은0회다. 공개 DFT에 대한 새 모델 평가와 정적 원자 이완은 실제 실행했다.','small')
    r.p('단위: 원자료의 길이 표기 A는 Å(ångström), 즉10<super>-10</super>m를 뜻한다. eV/A는 원자 힘, eV/A²는 에너지 곡률의 단위다.','small')
    if len(args.revision)==40 and all(c in '0123456789abcdef' for c in args.revision):
        revision_link=f'<link href="https://github.com/donghyeok8649daniel/al-fatigue-probability-theory/tree/{args.revision}" color="#087e8b">{args.revision[:12]}</link>'
    else:
        revision_link=escape(args.revision)
    r.p(f'<b>코드·원자료:</b> silicon-wafer-research / {revision_link}<br/>기록 시각 {escape(args.as_of)} KST. PDF와 사용 파일의 SHA256은 함께 생성한 manifest에 기록했다.','small')

    r.page('01. 균열 개시를 어떻게 구분하는가','과거 성장 연구를 개시 결과로 재명명하지 않는다.')
    r.table(['연구 대상','개시 연구에서의 위치'],[
        ['v2 crack seed 선단 전이, v3/v4 해당 끝점 열운동','기존 균열 주변 성장 관련 보조 자료. 0.703eV를 무균열 시편의 개시 장벽으로 쓰지 않는다.'],
        ['v6 K 경계와 약110MPa 환산','이미 균열이 있는 기하의 에너지 비교. 초기 시편의 최초 개시 응력이 아니다.'],
        ['v10 B/P 강체 계면 분리','조성·배치별 정적 분리 일 대조. 핵생성 경로·유한온도 개시 장벽이 아니다.'],
        ['v11 유한 Si 시편','내부 면 절단·공극·crack seed를 넣지 않았다. 외부 bare 표면과 rigid grip은 명시했다.']],[1.2,2.4])
    r.section('첫 균열의 정의에는 추가 검증이 필요하다')
    r.p('단일 결합이 길어짐, 원자 이웃 그래프의 변화, 응력의 최대값은 관측 진단이다. 이들 중 하나를 임의 임계값으로 잘라 첫 균열이라고 선언하지 않는다. 쉽게 되돌아오는 결합 개방과 지속되는 균열 상태를 구분해야 한다.')
    r.box('초기 상태 A와 균열 상태 B를 정의한 뒤, B에 처음 들어가는 확률을 계산한다. 그 B가 실제 균열 개시를 나타내는지와 그 확률을 계산하는 수치식이 맞는지는 별개의 검증이다.')
    r.p('표면, 계면, 도펀트 근처를 서로 다른 원인별 B로 나눌 수 있다. 원인별 집합은 겹치지 않아야 하며, 같은 사건을 두 번 흡수 질량으로 세지 않는다. 새로운 재료도 공통 확률 구조 안에서 다룬다.')
    r.p('추가 검토한 <link href="https://doi.org/10.1016/j.mtcomm.2022.105002" color="#087e8b">native-oxide Si nanowire 연구(2023)</link>의 p=(W0-Wfail)/W0는 파손 시 폭 감소 지표다. 이름에 probability가 있어도 ensemble의 첫 균열 흡수확률과 같지 않다. 원문 식과 그림을 확인해 구분했다.','small')

    r.page('02. 제공 논문이 측정한 것','Tsuchiya 외, Engineering Fracture Mechanics 163 (2016) 523-532')
    r.p('<link href="https://doi.org/10.1016/j.engfracmech.2015.08.029" color="#087e8b">Fracture behavior of single crystal silicon with thermal oxide layer</link>를 제공 PDF의 표·본문·그림으로 검토했다. 단조 인장 후 파단과 파면을 분석한 논문이다. 시간에 따른 첫 균열 관측이나 피로 전이율을 측정한 자료는 아니다.')
    r.table(['항목','원문의 조건/범위'],[
        ['시편','(100) SOI, 두께5μm·폭4μm, 길이120/600μm, 면내 &lt;110&gt; 인장.'],
        ['실험','26°C, 상대습도50%, stage 속도0.5μm/s. 도펀트 종·농도는 보고하지 않았다.'],
        ['열산화','1100°C, 5.5/22.5/90분. 예상 막두께50/100/200nm이며 개별 미세시편에서 직접 측정한 두께는 아니다.'],
        ['표면','가공 scallop 간격120nm·깊이30nm. 논문의 FEM은 scallop을 생략한다.'],
        ['FEM 입력','측정 평균 파단강도를 외력으로 사용. 산화막 E=70GPa, 참고문헌에서 가져온 -460MPa 잔류응력.']],[1,3])
    labels=['Bare','50nm','100nm','200nm','Oxide removed']
    r.table(['조건','측정 명목 파단강도(GPa)','논문 Si FEM 응력(GPa)'],[
        [name,f'{s["nominal_strength_GPa"]:.2f}',f'{s["FEM_Si_axial_GPa"]:.2f}'] for name,s in zip(labels,source['rows'])],[1.1,1.3,1.3])
    r.p('Si FEM 응력은 측정 강도를 하중으로 넣은 뒤 계산한 값이다. 이를 독립적인 파단강도 예측 검증으로 세면 순환 논리가 된다. 파면에서 추정한 기점도 시간 분해된 첫 개시 관측과 다르다.','small')

    r.page('03. 산화층의 단면과 응력 분담','논문 수치의 재현과 원문 기하 규약의 불일치를 분리했다.')
    r.picture('figures_stage1/thermal_oxide_audit.png',210)
    r.section('27조건을 따로 계산했다')
    r.p('외형치에서 면당 두께를 빼는 경우, 양면 전체 두께로 읽는 경우, 평면 산화에서 Si 소비비44%를 가정하는 경우를 구분했다. 잔류응력도 입력의 의미가 다른3가지로 나눴다. 기하를 바꿀 때 실험 외력은 고정했다.')
    r.p('원문 외형치와 면당 막두께 규약에서는 이론 탄성계수 차이가 최대0.103GPa, Si 응력은 원문 FEM값과 최대0.609% 차이다. 이는 응력 분담의 산술 대조다. 1D 힘 잔차는 최대1.42e-14GPa·μm²였다.')
    r.box('200nm 막을 면당 두께로 읽으면 Si core는3.71×4.71μm이다. 논문 산화막 제거 행은3.91×4.91μm여서 두 방향 모두0.20μm 차이가 남는다. 양면 합계 두께 해석은 제거 행과 맞지만 이론 탄성계수 재현이 나빠진다. 원문 규약은 미해결이다.')
    r.p('Eq.3의 반경은1100°C 조건의 산화 과정 가설이다. Eq.4에서 파단강도로 역산한19~49nm 결함 크기는 독립적으로 측정한 초기 결함 분포가 아니다. 같은 강도를 다시 예측하는 입력으로 재사용하지 않았다.','small')

    r.page('04. 외부 DFT로 재료 에너지를 감사','모델 구현의 일치와 물리적 적합성은 같은 판정이 아니다.')
    r.table(['원자료','평가 범위','실행/독립성'],[
        ['MTPu / QE PBE','1,159상태, 17,413원자. 전체 포함.','공개 학습 자료와 새 MACE 힘 비교. 새 DFT0.'],
        ['MLFF-SiOx / CP2K PBE','1,466개 중 N≤320인1,013개. 제외453개 ID보존.','H/O/Si를 그대로 평가. force error로 제외하지 않음.'],
        ['MACE-MP-0b3-medium','CPU float64, neutral model.','원래 MTP/GAP 재실행 아님. 사전학습 중복 미감사.']],[1,1.5,1.5])
    r.picture('figures_stage1/oxide_force_audit.png',220)
    r.p('QE SiO2 조성의 힘 성분 RMSE는0.15823eV/A, 다른 Si/O 조성은0.34272eV/A다. Si-only는 큰 오차 상태를 포함해1.27379eV/A다. 조성 정보만 있는 상태를 산화 계면으로 부르지 않았다.')
    r.p('독립 고정열 parser로 QE1159개 원본 배열과 raw 배열을 대조하고 프레임/그룹 지표 차이0을 확인했다. 유효한 left-handed cell을 처음에 거부한 parser 결함은 수정 후 전체 재실행했으며 실패본도 보존했다.','small')

    r.page('05. 산화 표면의 힘 오차','CP2K 자료1,013상태와179,857원자의 실제 새 모델 평가')
    r.picture('surface_environment_summary_v2/surface_force_groups.png',235)
    r.table(['그룹','상태/원자 관측수','힘 성분 RMSE(eV/A)'],[
        [g['group'],f'{g["frames"]} / {g["atoms"]:,}',f'{g["force_component_RMSE_eV_A"]:.6f}']
        for g in surface['groups'] if not g['group'].startswith('small_')],[1.3,1,1])
    r.p('9개 표준 MACE 대조에서 force-only 경로와 에너지·힘 차이는0이었다. 독립 shlex/고정열 parser 대조의 raw 배열 차이0, 프레임 지표 차이 최대3.35e-14, 그룹 지표 차이0을 확인했다.')
    r.box('vacuum1은 주기 cell에8A를 넘는 빈 방향이 하나 있다는 기하 진단이다. 개별 원자의 계면 소속·전하·균열 여부를 정한 라벨이 아니다. 작은 시스템의 큰 오차도 오른쪽 그림에 그대로 남겼다.')

    r.page('06. 어느 국소 환경을 먼저 개선할까','산소 이웃1~3개 Si의 큰 오차는 판정 거리를 바꿔도 유지된다.')
    r.picture('surface_environment_summary_v2/surface_chemical_environments.png',240)
    central=[g for g in env['oxidized_surface_chemical_classes'] if g['radius_set']==1]
    r.table(['기하 분류','원자비중','오차제곱합 비중','힘 RMSE(eV/A)'],[
        [g['category'],f'{100*g["atom_fraction"]:.2f}%',f'{100*g["squared_error_fraction"]:.2f}%',f'{g["force_component_RMSE_eV_A"]:.4f}'] for g in central],[1.3,.8,1,1])
    r.p('산소 이웃1~3개 Si의 RMSE는 Si-O거리1.8/2.0/2.2A에서0.5735/0.5720/0.5700eV/A다. 기하 분류의 인구가 달라도 큰 오차가 유지된다. 산화수나 최초 균열 위치를 증명한 것은 아니다.')
    r.p('산화 Si/O/H군에서 프레임별 균일 힘오차 성분은 총 오차제곱합의0.0006553%다. 원자료의 nonzero 총힘은 전체 모델 오차를 설명하지 못한다. 원래 힘과 RMSE를 보존했다.','small')

    r.page('07. 힘뿐 아니라 에너지 차이도 비교','같은 H/O/Si 원자수끼리 비교해 임의의 조성별 에너지 상수를 상쇄했다.')
    r.box('error(i,j) = [(E_ML,j - E_ML,i) - (E_DFT,j - E_DFT,i)] / N<br/>에너지 offset이나 배율을 적합하지 않았다. 상태쌍은 독립 표본이나 연속 핵생성 경로가 아니다.')
    r.table(['관찰 그룹','상태 수','동일조성 쌍','차이 RMSE(eV/atom)'],[
        [g['group'],g['frames'],g['pairs'],f'{g["pair_energy_error_RMSE_eV_atom"]:.6f}'] for g in energy['groups'] if not g['group'].startswith('small_')],[1.3,.6,.8,1.1])
    r.section('상태의 에너지 순서가 바뀌는 경우도 있다')
    r.p('산화 Si/O/H군에서 DFT차이0이 아닌25,954쌍 중2,107쌍은 모델 에너지 차이의 부호가 반대였다. DFT차이의 절댓값이1e-4eV/atom보다 큰 쌍만 표기하면2,056/25,863이다. 이는 저장 스냅샷의 순서 진단이며 파손 확률이나 독립 오분류율이 아니다.')
    r.p('기존 첫 동조성 anchor 기반 상대에너지와 raw 재계산의 차이는0이다. 전체 쌍 오차제곱합과 분산 항등식의 상대차이도 최대1.62e-13으로 확인했다. 작은 시스템의 큰 에너지 오차도 결과 표에 보존했다.')
    r.box('평균 에너지 오차가 작아 보여도 개시 장벽의 정확성을 보장하지 않는다. 실제 무균열 시작점과 후보 전이상태·끝점을 같은 전자구조 규약으로 비교하는 경로 검증이 필요하다.')

    r.page('08. 미리 균열을 넣지 않은 시편','360Si 원자, 자유원자216개, 양쪽 rigid grip144개. bare 외부 표면.')
    r.p('인장 방향은[110], 명목 단면은190.4845A², 초기 grip 중심 간격은30.9486A다. 외부 자유표면을 만들었지만 내부 균열·공극·잘린 면을 넣지 않았다. 산화막이나 도펀트가 들어간 실험 미세시편 전체는 아니다.')
    snapshot=args.snapshot or 'prism_snapshot_audit_0406'
    r.picture(snapshot+'/prism_geometry.png',280)
    r.p('그림은 독립 snapshot 감사가 지정한 마지막 수렴 상태다. 보라색은 grip, 주황색은 초기 격자 대비 새 거리이웃 쌍이다. 결합 재구성을 곧바로 균열로 세지 않았다.','small')
    r.table(['변위제어 strain','명목응력(GPa)','자유원자 최대 잔여힘(eV/A)','수렴'],[
        [f'{100*s["strain"]:.1f}%',f'{s["nominal_stress_GPa"]:.4f}',f'{s["free_force_max_eV_A"]:.3e}','예' if s['converged'] else '아니오'] for s in displacement[-6:]],[1,1,1.7,.5])
    r.p(f'표는 최근 최대6개 저장 상태다. 저장 {len(displacement)}개 중 수렴 {len(converged)}개. 3시간 예산으로10% 상태에서 중단됐다. 10%의 큰 잔여힘과 응력 하락은 미수렴 기록이며 파손강도 또는 개시 판정으로 쓰지 않는다.','small')
    continuation_path=root/'prism_10pct_continuation/summary.json'
    if continuation_path.exists():
        continuation=r.data('prism_10pct_continuation/summary.json')
        r.p(f'별도10% 후속 이완: 수렴 {"예" if continuation["converged"] else "아니오"}, 최대 자유힘{continuation["free_force_max_eV_A"]:.3e}eV/A, 응력{continuation["nominal_stress_GPa"]:.4f}GPa. 위 원래 중단 자료를 덮어쓰지 않고 같은 grip 조건에서 이어갔다.','small')
    elif (root/'prism_10pct_continuation/protocol.json').exists():
        r.data('prism_10pct_continuation/protocol.json')
        r.p('별도10% 후속 이완이 진행 중이다. 저장된 full-precision 좌표와 고정 grip을 그대로 이어받았으며3000초 한도다. 아직 수렴 판정은 없다.','small')

    r.page('보충. 인장 중 원자 재배열','전체 자유원자의 균일 변형과 이동을 제거한 기하 잔차')
    rearrangement=r.data(args.rearrangement+'/summary.json')
    r.picture(args.rearrangement+'/prism_nonaffine.png',260)
    r.table(['인장 구간','잔차 RMS(A)','최대 잔차(A)','grip 이웃의 제곱합 비중'],[
        [f'{100*s["strain_from"]:.0f}→{100*s["strain_to"]:.0f}%',f'{s["free_nonaffine_RMS_A"]:.5f}',f'{s["free_nonaffine_max_A"]:.5f}',f'{100*s["fraction_squared_residual_at_grip_neighbors"]:.2f}%'] for s in rearrangement['intervals']],[1,1,1,1.5])
    r.p('자유원자216개의3D affine 변환과 translation을 최소제곱으로 맞춘 뒤 남는 위치 차이다.6→8% 구간의 RMS는 앞선4→6%보다6.45배 크다. grip에 직접 이웃한36개보다 나머지180개에서 RMS가 크다.')
    r.box('큰 구조 재배열은 확인했지만, 이 값은 소성변형률이나 에너지 국소화 또는 균열 detector가 아니다. grip 이웃에만 집중되지 않는다고 해서 경계조건의 영향을 배제한 것은 아니다.')
    r.p(f'강체 회전·이동에 대한 잔차 크기 불변성 오차는 최대{rearrangement["maximum_rigid_frame_invariance_error_A"]:.2e}A다. 수렴한5상태만 사용했고 새 potential/DFT/MD 평가는0회다.','small')

    r.page('보충. 10%에서의 구조 변화','후속 이완은 수렴했지만 최초 균열 판정은 아직 하지 않았다.')
    r.picture('prism_10pct_replay/continuation_geometry.png',230)
    r.table(['관측량','후속 이완의 결과'],[
        ['추가 원자 이완',f'{continuation["force_calls"]} force calls / {continuation["elapsed_s"]:.2f}s 경과시간'],
        ['명목응력 / 최대 자유힘',f'{continuation["nominal_stress_GPa"]:.6f}GPa / {continuation["free_force_max_eV_A"]:.3e}eV/A'],
        ['원시 관측량 / grip 재검증',f'{continuation_audit["maximum_observable_replay_difference"]:.1e} / {continuation_audit["fixed_grip_difference_A"]:.1e}A'],
        ['8% 대비 affine 제거 잔차',f'RMS {continuation_audit["previous_to_current_nonaffine_RMS_A"]:.5f}A / 최대 {continuation_audit["previous_to_current_nonaffine_max_A"]:.5f}A']],[1.3,2])
    graph=continuation_audit['graphs'][0]
    r.p(f'2.8A 거리기준에서 초기 격자 이웃쌍{graph["lost_from_initial_lattice"]}개가 사라지고{graph["formed_from_initial_lattice"]}개가 새로 생겼다. 네 거리규약 모두 전체 원자 연결은 유지된다. 이는 큰 재배열의 진단이며 최초 균열의 인증은 아니다.')
    r.p(f'같은 고정10% 경계에서 optimizer가 에너지를{-continuation_audit["energy_change_at_fixed_displacement_eV"]:.5f}eV 낮췄다. 이 구간 grip의 추가 일은0이다. 정적 최적화에는 물리 시간이 없으므로 이 값을 물리 소산이나 전이 장벽으로 부르지 않는다.','small')
    if (root/'prism_unload_replay/summary.json').exists():
        unload=r.data('prism_unload_replay/summary.json')
        r.box(f'8%로 되돌린 이완의 수렴: {"예" if unload["unloading_converged"] else "아니오"}, 자유힘 최대{unload["unloading_free_force_max_eV_A"]:.2e}eV/A. 이전8%와 같은 grip에서 에너지 차이{unload["energy_difference_at_same_grips_eV"]:.5f}eV, 응력{unload["loading_stress_GPa"]:.4f}→{unload["unloading_stress_GPa"]:.4f}GPa다. 동일 Si 원자 순서에 무관한 최소 matching RMS는{unload["minimum_matching_RMS_difference_A"]:.5f}A다. 정적 하중 이력 대조이며 실제 피로주기 또는 물리 소성 검증은 아니다.')
    else:
        r.box('이10% 구조를8%로 되돌리는 별도 정적 이완이 진행 중이다. 이전8%와 같은 grip에서 에너지·힘·원자 순서에 무관한 구조 차이를 대조할 예정이다. 아직 복원 여부의 결론은 없다.')

    r.page('보충. 같은8%로 되돌린 결과','처음8%와 경계 좌표가 같아도 정적 이완 구조와 반력은 달랐다.')
    unload=r.data('prism_unload_replay/summary.json')
    unload_run=r.data('prism_unload_10_to_8/summary.json')
    r.data('prism_unload_figure_v2/summary.json')
    r.picture('prism_unload_figure_v2/static_return.png',250)
    r.table(['같은8% 경계에서 비교','최초 인장','10%를 거친 되돌림'],[
        ['명목 축응력(GPa)',f'{unload["loading_stress_GPa"]:.6f}',f'{unload["unloading_stress_GPa"]:.6f}'],
        ['에너지(eV)',f'{unload["loading_energy_eV"]:.6f}',f'{unload["unloading_energy_eV"]:.6f}'],
        ['자유원자 힘 최대(eV/A)',f'{unload["loading_free_force_max_eV_A"]:.3e}',f'{unload["unloading_free_force_max_eV_A"]:.3e}']],[1.4,1,1])
    r.p(f'새 정적 이완{unload_run["force_calls"]}force calls/{unload_run["elapsed_s"]:.2f}s, {unload_run["attempts"][0]["steps"]}steps로 수렴했다. 같은 grip 좌표의 차이와 원시 관측량 재계산 차이는 모두0이다. 동일 Si 원자의 순서 교환까지 허용한 최소 거리 RMS는{unload["minimum_matching_RMS_difference_A"]:.5f}A이며, 최적 대응에서 바뀐 원자 번호는{unload["minimum_matching_permuted_atoms"]}개다.')
    r.p('초기8% 상태와2.8A 거리 이웃을 비교하면41쌍이 사라지고70쌍이 생겼다. 네 거리규약에서 전체 원자 연결은 유지된다. 그림의 선은 저장 끝점의 순서를 표시하며 연속적인 물리 하중 궤적을 측정한 것은 아니다.','small')
    r.box('선택한 작은 bare 시편과 정적 최적화 절차에서 이전8% 구조로 돌아오지 않았다. 에너지 차이-9.71922eV는 같은 경계의 두 정적 상태 차이다. 실제 사이클의 소산 일·물리 소성·균열 개시 확률로 해석하지 않는다. 두 상태의 국소 안정성 및 관찰된 재배열의 물리적 종류는 추가 검증이 필요하다.')

    r.page('09. 변위0과 축력0은 다르다','0/100MPa 힘제어는 U - F·Delta를 자유원자와 grip 연장에 대해 최소화했다.')
    r.picture(args.force_replay+'/force_control_response.png',240)
    r.table(['목표 축응력','계산 축응력(MPa)','축력 오차(MPa)','전체 지지잔차/면적(MPa)'],[
        [f'{1000*s["target_nominal_GPa"]:.0f}MPa',f'{1000*s["actual_nominal_GPa"]:.6f}',f'{s["axial_target_residual_MPa"]:.6f}',f'{s["support_resultant_equivalent_MPa"]:.6f}'] for s in force['records']],[.9,1.1,1.1,1.4])
    comparison=force['comparison']
    r.p(f'초기 길이를 고정한 상태는 표면 이완 후3.184GPa의 인장 반력을 가진다. 별도로 축력0을 풀면 grip 간격이 초기값보다{-100*comparison["zero_load_contraction_relative_to_bulk_grip_length"]:.5f}% 줄어든다. 100MPa에서는 이 축력0 길이 대비{100*comparison["strain_relative_to_zero_axial_load"]:.7f}% 늘어난다.')
    r.p(f'두 점의 응력/변형률 secant는{comparison["secant_stress_strain_ratio_GPa"]:.3f}GPa다. 매우 작은 bare 시편과 rigid grip 조건의 응답이며 bulk 영률이나 실험 웨이퍼의 검증값이 아니다. 원시 좌표·힘·gradient·enthalpy를 독립 재계산했다.')
    if args.force_replay=='force_control_replay_tight':
        tolerance=r.data('force_tolerance_comparison/summary.json')
        r.p('힘 허용오차를5e-4에서1e-5eV/A로50배 엄격하게 한 대조다. 전체 지지잔차/면적은 무하중5.049→0.071MPa,100MPa 하중0.823→0.064MPa로 줄었다. 두 점의 secant도112.438→112.163GPa로 바뀌었다. 원래 결과는 보존했다.','small')
    r.box('축력 목표 오차가 작아도 모든 지지 반력이 정확히 평형인 것은 아니다. 표의 마지막 열은 유한 수렴 오차다. 이를0으로 숨기거나 실제 물리 외력으로 재해석하지 않았다.')

    r.page('10. 국소 안정성과 경계 일의 검증','힘 수렴, 국소 곡률, 전역 안정성, 물리 개시는 서로 다른 단계다.')
    audit_path=root/args.curvature_audit/'summary.json'
    if audit_path.exists():
        audit=r.data(args.curvature_audit+'/summary.json')
        r.table(['항목','실제 계산 상태'],[['Lanczos 수렴','수렴' if audit['lanczos_converged'] else '미수렴'],
            ['중단/미수렴 기록',escape(str(audit['interruption'])) if audit['interruption'] else '없음'],
            ['검사된 mode 수',len(audit['modes'])],
            ['U gradient 최대(eV/A)',f'{audit["U_reduced_gradient_max_eV_A"]:.6e}'],
            ['U-F·Delta gradient 최대(eV/A)',f'{audit["enthalpy_reduced_gradient_max_eV_A"]:.6e}']],[1.5,1])
        if audit['modes']:
            r.table(['Mode','Ritz 곡률(eV/A²)','반간격 Rayleigh','잔차(eV/A²)'],[
                [m['mode'],f'{m["ritz_curvature_eV_A2"]:.6e}',f'{m["half_step_rayleigh_eV_A2"]:.6e}',f'{m["residual_norm_eV_A2"]:.3e}'] for m in audit['modes']],[.5,1.2,1.2,1])
    else:
        r.box('작성 시점에는 국소 곡률과 힘제어 후처리의 최종 결과가 아직 없다. 아래 내용은 실행 중 검사의 정의이며 통과 판정이 아니다.')
    r.section('계산하는 허용 변위 공간')
    r.p('고정 변위에서는 자유원자648좌표를 검사한다. 힘제어에서는 grip 연장1좌표까지 포함한649좌표를 사용한다. 외력의 일F·Delta가 선형이므로 에너지의 둘째 미분은 같지만, 평형을 판단하는 첫 미분에서는 F를 빼야 한다.')
    r.p('첫100MPa/Krylov16 검사는75분/519force calls에서 미수렴 종료됐고 mode0개다. 이를 보존한 뒤, 더 엄격하게 이완한100MPa 상태와 Krylov64로 별도 대조했다. 원자 상태와 Krylov 크기가 함께 달라져 단일변수 수렴 검사는 아니다.','small')
    if audit_path.exists() and audit.get('krylov_dimension')==64:
        check_note=('반간격 Hessian과 ±에너지 차분은 원자료의 mode별 검사에 기록했다.'
                    if audit['independently_checked_modes'] else
                    '독립 검사 mode가 없어 반간격·에너지 차분에 근거한 안정성 판정은 하지 못했다.')
        r.p(f'표는 정밀100MPa/Krylov64 대조다. 실제{audit["force_calls"]}force calls/{audit["elapsed_s"]:.2f}s, 저장mode{audit["saved_eigenmodes"]}개/독립 검사{audit["independently_checked_modes"]}개. {check_note} 기울기 raw 재검증 차이{audit["maximum_force_state_gradient_replay_difference_eV_A"]:.2e}eV/A다.','small')
    else:
        r.p('위 표는 첫 검사의 기록이다. 정밀100MPa/Krylov64 후처리 결과는 아직 이 보고서에 포함되지 않았다. 첫 검사에서는 수렴mode가 없어 반간격·독립 에너지 검사를 실행하지 못했다.','small')
    r.box('안정하게 재현되는 음의 곡률은 국소 불안정성의 근거가 될 수 있다. 제한된 Lanczos에서 나온 양의 Ritz값만으로 모든 방향의 양성이나 전역 최소를 증명할 수는 없다. 원자 질량에서 물리 clock을 만들어 넣지 않는다.')

    r.page('11. 도핑과 공통 첫 통과 확률','도펀트 원자수, 활성전하, 초기 화학 배치의 통계를 구분한다.')
    r.p(f'현재360자리 시편의 명목 체적은{occupancy["nominal_volume_cm3"]:.5e}cm³다. 도펀트1개는{occupancy["one_dopant_per_prism_concentration_cm3"]:.5e}cm<super>-3</super>에 해당한다. 따라서 낮은 농도를 모든 작은 상자에 도펀트1개씩 넣는 방식으로 표현할 수 없다.')
    r.box('균일하고 독립적인 치환을 가정하면 C=1e18cm<super>-3</super>에서 이 시편에 도펀트가 적어도1개 있을 확률은0.7342%다. 실제 편석·활성전하·상관 배치를 측정한 값이나 균열 개시확률은 아니다.')
    r.section('고정된 화학 배치 D를 먼저 조건으로 둔다')
    r.p('S_mix(t) = sum_D w_D S_D(t)<br/>생존 시편에서의 배치 가중치 = w_D S_D(t) / S_mix(t)<br/>평균 에너지나 평균 전이율을 한 번 대입하는 것과 일반적으로 다르다. 각 배치에서 같은 확률 구조를 사용하며 물질별 별개 확률 이론을 만들 필요는 없다.')
    r.section('실제로 구현·검증한 공통 수학')
    r.p('dP/dt = L P의 transient block과 원인별 흡수 collector를 분리했다. 생존 질량과 원인별 누적 첫 통과 확률을 함께 보존한다. committor는 A로 되돌아가기 전에 B에 도달할 확률이며, MFPT는 정의된 B까지의 평균 첫 도달시간이다.')
    r.p('겹치는 원인집합, 흡수되지 않는 닫힌 상태군, 유효하지 않은 generator를 거부한다. 해상도가 부족한 희귀 흡수 문제를 정칙화나 clipping으로 유한한 물리 보정값으로 만들지 않는다.')
    r.p('이번 Si/O/H 원자료에는 B/P가 없다. 위 산화 표면 대조를 도핑 계면 검증으로 부르지 않는다. v10 B/P 중성 강체분리 결과도 실제 활성전하별 개시 장벽을 대신하지 않는다.','small')

    r.page('12. 검증, 실패, 남은 연구','수치 구현의 완료와 재료·물리 시간의 미완료를 함께 기록한다.')
    r.picture('figures_stage1/first_passage_validation.png',200)
    r.p(f'관련 구현 {validation["tests_passed"]}PASS ({validation["test_elapsed_seconds"]:.2f}s). 합성 예제 최대 질량 잔차{math["maximum_mass_residual"]:.2e}, 확산 격자 오차비{math["diffusion_error_ratios"][0]:.2f}/{math["diffusion_error_ratios"][1]:.2f}. 실제 Si의 전이율·초·Hz를 검증한 값은 아니다.','small')
    r.table(['보존한 실패/한계','처리'],[
        ['초기 LBFGS 정체','미수렴 checkpoint 보존 후 line-search 재시작.'],
        ['extxyz 자릿수 차이','입력 저장오차 확인 후 prescribed grip 좌표를 정확 복구. 힘 허용오차를 늘린 수정이 아님.'],
        ['left-handed cell parser','abs(det)>0 규약 수정, 동등 기저 대조 및1159상태 새 전체 실행.'],
        ['변위제어3시간 예산 종료','0~8%5상태 수렴.10% 미수렴 구조·힘·에너지는 그대로 보존.'],
        ['미해결 물리','개시 basin·핵생성 경로·DFT 장벽·유한온도 generator·이동도·physical time.']],[1,2.5])
    manifest_status='있음' if (root/'artifact_manifest.json').exists() else '작성 중'
    r.p(f'<b>재현 위치:</b> results/silicon_initiation_v11의 REPRODUCE.md, WORKING_STATUS.md, validation.json, 각 raw/protocol/summary. 최종 source/artifact manifest는{manifest_status}. 이 PDF 옆 manifest는 사용한 파일의 실제 SHA256과 페이지 경계를 기록한다.','small')
    r.p('<b>주요 원자료:</b> <link href="https://doi.org/10.1016/j.engfracmech.2015.08.029" color="#087e8b">Tsuchiya 외(2016)</link>; <link href="https://doi.org/10.1038/s41524-024-01390-8" color="#087e8b">Zongo 외(2024), MTPu</link>; <link href="https://doi.org/10.1063/5.0220091" color="#087e8b">Cvitkovich 외(2024), MLFF-SiOx</link>. 세부 원문/초록 확인 범위는 INITIATION_SOURCE_SCOPE.md에 기록했다.','small')
    r.p('다음 핵심은 무균열 초기 상태와 지속되는 균열 상태의 경로를 정의하고 해당 Si/O/도펀트 환경의 에너지·힘을 검증하는 것이다. 현 단계에서 GPa를 MPa로 임의 축소하거나 정적 장벽에 임의 clock을 붙여 수명을 제시하지 않는다.','small')
    r.close();print(json.dumps(dict(output=args.output.name,pages=r.n,source_files=len(r.used),label=args.label)))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--results',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--font',type=Path,required=True);p.add_argument('--bold-font',type=Path,required=True)
    p.add_argument('--as-of',required=True);p.add_argument('--label',default='진행 중 / 중간 기록')
    p.add_argument('--revision',default='uncommitted v11 research');p.add_argument('--snapshot')
    p.add_argument('--force-replay',default='force_control_replay')
    p.add_argument('--rearrangement',default='prism_rearrangement_0610')
    p.add_argument('--curvature-audit',default='curvature_force_ensemble_audit')
    main(p.parse_args())
