import pytest

from .report_core_validation_ledger import read_report


def test_executed_report_counts_do_not_turn_skips_into_passes(tmp_path):
    path=tmp_path/'actual.xml'
    path.write_text('<testsuites><testsuite tests="7" failures="0" errors="0" skipped="2" time="1.25"/></testsuites>',encoding='utf8')
    result=read_report(path)
    assert result['passed']==5 and result['skipped']==2
    assert result['elapsed_seconds']==1.25
    assert len(result['sha256'])==64


def test_failed_or_empty_reports_cannot_pass(tmp_path):
    path=tmp_path/'actual.xml'
    for content in ('<testsuite tests="2" failures="1"/>',
                    '<testsuite tests="0"/>','<testsuites/>'):
        path.write_text(content,encoding='utf8')
        with pytest.raises(ValueError):read_report(path)
