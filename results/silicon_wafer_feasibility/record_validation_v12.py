"""Run the relevant implementation checks once and bind their exact source bytes."""
from __future__ import annotations
import argparse,hashlib,json,re,subprocess,sys
from datetime import datetime,timezone
from pathlib import Path

MODULES=('silicon_charge_dynamics','silicon_initiation_probability','silicon_initiation_research',
         'silicon_structure_diagnostics_v12','silicon_oxide_delta_v12','silicon_oxide_angular_delta_v12')


def main(args):
    repo=Path(__file__).resolve().parents[2]
    if args.output.exists():raise ValueError('do not overwrite completed validation evidence')
    tests=['solver_v1/test_'+name+'.py' for name in MODULES]
    sources=['solver_v1/'+name+'.py' for name in MODULES]+tests+[
        'solver_v1/silicon_oxide_research.py','solver_v1/probability_pde_2d.py']
    before={p:hashlib.sha256((repo/p).read_bytes()).hexdigest() for p in sources}
    utc=datetime.now(timezone.utc).isoformat()
    result=subprocess.run([sys.executable,'-m','pytest',*tests,'-q'],cwd=repo,capture_output=True,text=True)
    log=result.stdout+'\n'+result.stderr
    # Preserve warning text while excluding personal absolute runtime paths
    # from the public research bundle. This is not a warning suppression.
    for absolute,label in sorted(((str(repo),'<REPO>'),(str(Path(sys.prefix)),'<PYTHON_ENV>')),
                                 key=lambda pair:len(pair[0]),reverse=True):
        log=log.replace(absolute,label).replace(absolute.replace('\\','/'),label)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.with_name('implementation_tests.txt').write_text(log,encoding='utf-8')
    if result.returncode!=0:raise RuntimeError('implementation tests failed; exact output preserved')
    match=re.search(r'(\d+) passed(?:, \d+ warnings?)? in ([\d.]+)s',log)
    if not match:raise ValueError('actual pytest summary unavailable; no synthetic pass count')
    for path,sha in before.items():
        if hashlib.sha256((repo/path).read_bytes()).hexdigest()!=sha:
            raise ValueError('source changed while tests were running: '+path)
    evidence=dict(tests_passed=int(match.group(1)),test_elapsed_seconds=float(match.group(2)),
        command='python -m pytest '+' '.join(tests)+' -q',log_summary=match.group(0),
        log_path_redaction='repository and Python environment absolute paths replaced by labels; warning text retained',
        started_UTC=utc,finished_UTC=datetime.now(timezone.utc).isoformat(),
        source_files=[dict(path=p,sha256=sha) for p,sha in sorted(before.items())],
        scope='probability conservation/first passage/charge dynamics; prism work and geometry; conservative oxide force derivatives and invariances',
        new_atomic_calculations=0,material_approved=False,kinetic_clock_calibrated=False)
    args.output.write_text(json.dumps(evidence,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(evidence,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args())
