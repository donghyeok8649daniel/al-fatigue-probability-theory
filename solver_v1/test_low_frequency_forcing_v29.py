import json
import numpy as np
import pytest
from .run_low_frequency_forcing_v29 import make_cases,validate_case,analyze


def test_predeclared_balanced_protocol_and_cycles():
    cases=make_cases([3.,1.5])
    assert len(cases)==12
    assert sum(c['duration_ps'] for c in cases)==10000
    for c in cases:
        assert c['frequency_per_ps']==.2
        assert (c['duration_ps']-100)*.2==round((c['duration_ps']-100)*.2)
        other=[x for x in cases if x['axis']==c['axis'] and x['fraction']==c['fraction']
            and x['dt_ps']==c['dt_ps'] and x['sign']==-c['sign']]
        assert len(other)==1
        assert other[0]['force_eV_A']==-c['force_eV_A']
        if c['fraction']==2: assert c['duration_ps']==500


@pytest.mark.parametrize('bad', [[0,1],[np.nan,1],[-1,2],[1]])
def test_bad_force_scale_refused(bad):
    with pytest.raises(ValueError):make_cases(bad)


def test_completed_case_binding_and_incomplete_rejected(tmp_path):
    case=make_cases([3.,1.5])[0]
    protocol=dict(restart_sha256='bound',frame_ps=.025)
    meta=dict(completed=True,restart_sha256='bound',duration_ps=case['duration_ps'],
        dt_ps=case['dt_ps'],frame_ps=.025,ensemble='nve',conjugate_drive=dict(
            axis=case['axis'],force_eV_A=case['force_eV_A'],frequency_per_ps=.2))
    path=tmp_path/'summary.json';path.write_text(json.dumps(meta))
    assert validate_case(case,protocol,tmp_path)['completed']
    for key,value in [('completed',False),('dt_ps',.005),('duration_ps',500),('restart_sha256','wrong')]:
        changed=dict(meta,**{key:value});path.write_text(json.dumps(changed))
        with pytest.raises(ValueError):validate_case(case,protocol,tmp_path)


def test_analysis_does_not_overwrite_existing_output(tmp_path):
    out=tmp_path/'results';out.mkdir()
    marker=out/'user.txt';marker.write_text('preserve')
    with pytest.raises(FileExistsError):analyze(tmp_path,tmp_path,out)
    assert marker.read_text()=='preserve'
