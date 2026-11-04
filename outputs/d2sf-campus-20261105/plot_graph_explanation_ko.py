"""Plot saved numerical results with Korean labels; no solver or fitting runs."""
from pathlib import Path
from fractions import Fraction
import json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.font_manager import FontProperties

ROOT=Path(__file__).resolve().parents[2]
SRC=ROOT/'outputs/sg-coarse-memory-20261104';OUT=Path(__file__).parent
r=json.loads((SRC/'CORRELATION_MEMORY_RESULT.json').read_text(encoding='utf-8'))
w=json.loads((SRC/'CORRELATION_WITNESS_EXACT.json').read_text(encoding='utf-8'))
prep=json.loads((SRC/'CORRELATION_PREPARATION_EXACT.json').read_text(encoding='utf-8'))
plt.rcParams.update({'font.family':'Malgun Gothic','axes.unicode_minus':False,'font.size':10,'pdf.fonttype':42,'ps.fonttype':42})
eight=next(x for x in r['comparisons'] if x['cycles']==8)
base=eight['summaries']['marginals_energy'];exact=eight['summaries']['plus_one_cycle_event_memory']
noise=next(x for x in r['noisy_probe_comparisons'] if x['cycles']==8 and x['probe_absolute_uncertainty']==1e-4)
fig,axs=plt.subplots(1,2,figsize=(11.4,4.4),gridspec_kw={'width_ratios':[.9,1.3]})
ax=axs[0]
vals=np.array([float(w[k]['decimal_lower']) for k in ['corrected_exact_cdf_plus','corrected_exact_cdf_minus']])*100
ax.bar(['준비집단 +','준비집단 -'],vals,color=['#555555','#aaaaaa'],width=.48)
ax.set_ylim(0,.27);ax.set_ylabel('한 주기 첫 열림 확률 (%)');ax.set_title('(가) 같은 요약값의 두 집단',pad=13)
for i,v in enumerate(vals):ax.text(i,v+.007,f'{v:.5f}%',ha='center',fontsize=11)
ax.text(.5,.94,'두 좌표의 전체 주변분포·평균 에너지는 같음',transform=ax.transAxes,ha='center',fontsize=8.5)
ax=axs[1]
for y,record,color in [(2,base,'#555555'),(1,noise,'#888888'),(0,exact,'#1d6271')]:
 lo,hi=100*record['lower'],100*record['upper']; mid=(lo+hi)/2
 ax.errorbar(mid,y,xerr=[[mid-lo],[hi-mid]],fmt='o',color=color,elinewidth=5,capsize=6,markersize=4)
 ax.text(mid,y+.20,f'{lo:.5f} ~ {hi:.5f}%',ha='center',fontsize=9)
ax.set_yticks([2,1,0],['현재 요약만','한 주기 응답 추가\n(오차 ±0.0001)','한 주기 응답 추가\n(오차 없음)'])
ax.set_xlim(.52,.765);ax.set_ylim(-.55,2.6);ax.set_xlabel('8주기 첫 열림 확률의 가능한 범위 (%)')
ax.set_title('(나) 추가 정보에 따른 식별범위',pad=13)
for ax in axs:
 ax.spines[['top','right']].set_visible(False);ax.grid(axis='x' if ax is axs[1] else 'y',alpha=.18)
fig.subplots_adjust(left=.075,right=.98,bottom=.19,top=.84,wspace=.62)
fig.text(.5,.035,'저장된 유한 SG 결과의 재표시. 두 패널은 서로 다른 준비집단 비교이며 실제 피로수명 실험이 아니다.',ha='center',fontsize=9)
fig.savefig(OUT/'APPLICATION_GRAPH_KO.png',dpi=240)
plt.close(fig)

ids=np.asarray(prep['indices']);shape=(len(prep['a']),len(prep['s']))
arrays=[]
for key in ['plus','minus']:
 a=np.full(np.prod(shape),np.nan);a[ids]=[100*float(Fraction(int(n),int(d))) for n,d in prep[key]];arrays.append(a.reshape(shape))
fig,axs=plt.subplots(1,2,figsize=(9.4,3.5),sharex=True,sharey=True)
vmax=max(np.nanmax(a) for a in arrays)
cmap=plt.get_cmap('Greys').copy();cmap.set_bad('#c9ced4')
for ax,a,title in zip(axs,arrays,['준비집단 +','준비집단 -']):
 im=ax.imshow(np.ma.masked_invalid(a),origin='lower',aspect='auto',norm=LogNorm(vmin=min(np.nanmin(x) for x in arrays),vmax=vmax),cmap=cmap,interpolation='nearest')
 ax.set_title(title);ax.set_xlabel('미끄럼좌표 s의 격자 번호')
axs[0].set_ylabel('열림좌표 a의 격자 번호')
fig.subplots_adjust(left=.07,right=.84,bottom=.23,top=.85,wspace=.12)
cax=fig.add_axes([.87,.25,.018,.53]);cb=fig.colorbar(im,cax=cax,label='초기 확률질량 (%, 로그 색상)');cb.set_ticks([100/150,1,5,10,25],labels=['0.667','1','5','10','25'])
fig.text(.45,.035,'회색은 준비영역 밖. 영역 안의 모든 칸은 공통 성분 0.667% 이상을 가진다. 행합·열합·평균 에너지는 같다.',ha='center',fontsize=8.5)
fig.savefig(OUT/'PREPARATION_HEATMAPS_KO.png',dpi=220);plt.close(fig)

source_files=['CORRELATION_MEMORY_RESULT.json','CORRELATION_WITNESS_EXACT.json','CORRELATION_PREPARATION_EXACT.json','REACHABLE_LOAD_ORDER_RESULT.json']
manifest={'source_sha256':{n:hashlib.sha256((SRC/n).read_bytes()).hexdigest() for n in source_files},'new_scientific_runs':0,'graph_kind':'Saved results only; first passage probability and identification interval, not physical fatigue accuracy','panel_a_probabilities':(vals/100).tolist(),'panel_b_horizon':8,'panel_b_intervals':{'summary_only':base,'noisy_response':noise,'exact_response':exact},'panels_use_distinct_preparation_families':True}
(OUT/'GRAPH_DATA_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print('KOREAN_GRAPH_AND_PREPARATION_MAPS_WRITTEN')
