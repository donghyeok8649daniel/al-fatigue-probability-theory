"""Capture exit status and full verification output, redacting local paths."""
import argparse
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

p=argparse.ArgumentParser()
p.add_argument('--output',type=Path,required=True)
p.add_argument('--name',required=True)
p.add_argument('command',nargs=argparse.REMAINDER)
a=p.parse_args()
command=a.command[1:] if a.command[:1]==['--'] else a.command
a.output.mkdir(parents=True,exist_ok=True)
start=time.perf_counter()
r=subprocess.run([sys.executable,*command],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
elapsed=time.perf_counter()-start
raw=r.stdout.decode('utf-8',errors='replace')
prefixes=[(str(Path.cwd()),'<REPO>'),(str(Path.home()),'<USER_PROFILE>')]
for old,new in prefixes:
    raw=raw.replace(old,new).replace(old.replace('\\','/'),new)
(a.output/(a.name+'.log')).write_text(raw,encoding='utf-8')
report=dict(command=['python',*command],exit_code=r.returncode,elapsed_reported_s=elapsed,
            python=platform.python_version(),path_redactions=['repository root','user profile'],
            timing_note='reported perf_counter duration; not inferred CPU time')
(a.output/(a.name+'.json')).write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
print(raw[-2500:]);print(json.dumps(report))
raise SystemExit(r.returncode)
