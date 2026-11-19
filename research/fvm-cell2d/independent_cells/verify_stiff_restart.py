"""New2D stiff action and actual disk checkpoint restart checks."""
import json,os,time
from pathlib import Path
import numpy as np
import fvm2d as m


m.set_model('coherent_g4')
folder=m.ROOT/'outputs/fvm-cell2d/component_restart'
folder.mkdir(parents=True,exist_ok=True)
grid=m.Grid2D(np.linspace(.7,6.,32),np.linspace(-.08,.25,14))
p,_=m.initial_mass(grid);mass=np.stack([p,np.roll(p,2)]);F=np.zeros(2)
op=m.assemble(grid,4.5);start=time.perf_counter()
explicit,ef,ea=m.advance_uniform_batch(op,mass,F,.03)
implicit,imf,ima=m.advance_stiff_batch(op,mass,F,.03)
L1=float(abs(explicit-implicit).sum(axis=1).max()+abs(ef-imf).max())
CDF=float(abs(ef-imf).max())
nextop=m.assemble(grid,1.)
uninterrupted,uf,_=m.advance_batch(nextop,explicit,ef,.01)
checkpoint=folder/'checkpoint.npz'
with checkpoint.open('wb') as file:np.savez_compressed(file,mass=explicit,F=ef,shape=grid.shape)
with np.load(checkpoint,allow_pickle=False) as saved:
    reloaded=np.asarray(saved['mass']);reloadedF=np.asarray(saved['F'])
resumed,rf,_=m.advance_batch(nextop,reloaded,reloadedF,.01)
restart_error=float(abs(resumed-uninterrupted).sum()+abs(rf-uf).sum())
report=dict(accepted=L1<2e-12 and CDF<2e-15 and restart_error==0.,
    active_dimension=2,model=m.model_metadata(),L1_BDF_vs_positive_uniformization=L1,
    CDF_BDF_vs_positive_uniformization=CDF,restart_difference=restart_error,
    checkpoint_bytes=checkpoint.stat().st_size,explicit_audit=ea,implicit_audit=ima,
    seconds=time.perf_counter()-start)
(folder.parent/'STIFF_AND_RESTART_CHECK.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2),flush=True)
assert report['accepted']
