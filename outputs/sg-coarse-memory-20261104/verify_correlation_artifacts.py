"""Re-open saved exact witness independently of the construction runner."""
from fractions import Fraction as F
from pathlib import Path
import json
import sys
sys.set_int_max_str_digits(50000)
p=Path(__file__).parent
j=json.loads((p/'CORRELATION_PREPARATION_EXACT.json').read_text())
plus=[F(*map(int,x)) for x in j['plus']]
minus=[F(*map(int,x)) for x in j['minus']]
energy=[F.from_float(float.fromhex(x)) for x in j['energy_binary64_hex']]
assert sum(plus)==sum(minus)==1 and min(plus+minus)>0
assert sum((x-y)*z for x,y,z in zip(plus,minus,energy))==0
ns=len(j['s']);ids=j['indices']
for a in range(len(j['a'])):
    assert sum((plus[k]-minus[k] for k,i in enumerate(ids) if i//ns==a),F(0))==0
for s in range(ns):
    assert sum((plus[k]-minus[k] for k,i in enumerate(ids) if i%ns==s),F(0))==0
w=json.loads((p/'CORRELATION_WITNESS_EXACT.json').read_text())
def value(k):return F(int(w[k]['numerator']),int(w[k]['denominator']))
assert value('corrected_exact_cdf_plus')-value('corrected_exact_cdf_minus')==value('exact_gap')>0
assert value('any_same_summary_predictor_max_error_lower_bound')*2==value('exact_gap')
c=json.loads((p/'CORRELATION_MEMORY_RESULT.json').read_text())
h=json.loads((p/'REACHABLE_LOAD_ORDER_RESULT.json').read_text())
assert c['checks_all_passed'] and h['checks_all_passed']
result={'status':'PASS_SAVED_EXACT_WITNESS_AND_COMPARISON_READBACK',
        'exact_positive_unit_mass':True,'exact_both_marginals_and_mean_energy':True,
        'exact_saved_cdf_gap_and_error_lower_bound':True,'LP_settings':16,
        'noisy_probe_settings':9,'load_protocols':len(h['cases']),
        'initial_readback_processing_failure':{'toolChunk':'2f3426','exit':1,
          'reason':'Default4300-digit integer conversion limit while reading generated5639-digit exact certificate',
          'resolution':'Bounded50000-digit local generated-certificate reader; mathematical inputs unchanged'},
        'independent_reviewer':False}
(p/'CORRELATION_ARTIFACT_VERIFY.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result))
