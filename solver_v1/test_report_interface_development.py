import pytest

from .report_interface_development import summarize_replay


def test_report_excludes_heldout_and_exact_constraints_from_loss():
    def row(name,role,r):return dict(observable=name,role=role,target='2',
        prediction=str(2+3*r),scale='3',normalized_residual=str(r))
    result=summarize_replay([row('bulk','exact',0.),row('force','fit',2.),
        row('direct_H','heldout',10.),row('bloch_q1','heldout',5.)])
    assert result['replayed_fit_loss']==4.
    assert result['interface_heldout_normalized_rms']==10.
    assert result['bloch_heldout_normalized_rms']==5.
    assert result['exact_count']==1


def test_report_refuses_inconsistent_scale_or_residual():
    row=dict(observable='x',role='fit',target='2',prediction='3',scale='2',normalized_residual='1')
    with pytest.raises(ValueError,match='residual'):summarize_replay([row])
    row['scale']='0'
    with pytest.raises(ValueError,match='positive'):summarize_replay([row])
