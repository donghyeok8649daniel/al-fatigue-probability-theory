from pathlib import Path
from datetime import datetime,timezone
import hashlib,json
import numpy as np

original=Path(__file__).resolve().parents[2]
si=original.parent/'aft-silicon-wafer'
root=si/'results/silicon_initiation_v12'
stage=Path(__file__).parent/'stage/results/silicon_initiation_v12'
stage.mkdir(parents=True,exist_ok=True)
events=json.loads((Path(__file__).parent/'power_events_gap.json').read_text(encoding='utf-8-sig'))
wake=next(e for e in events if e['provider']=='Microsoft-Windows-Power-Troubleshooter' and e['id']==1)
sleep=datetime.fromisoformat(wake['fields']['SleepTime'].replace('Z','+00:00'))
awaken=datetime.fromisoformat(wake['fields']['WakeTime'].replace('Z','+00:00'))
start=datetime.fromisoformat('2026-09-24T17:40:54+00:00')
deadline=datetime.fromisoformat('2026-09-25T01:40:54+00:00')
def dump(path,value):path.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
interruption=dict(status='interrupted_by_battery_hibernation; reporting after resume',
    requested_start_UTC=start.isoformat(),requested_deadline_UTC=deadline.isoformat(),
    requested_wall_hours=8,continuous_eight_hour_run_completed=False,
    Windows_sleep_UTC=wake['fields']['SleepTime'],Windows_wake_UTC=wake['fields']['WakeTime'],
    Windows_kernel_power_42_reason='Battery',Windows_resume_507_reason='Resume from Hibernate',
    reported_hibernation_wall_seconds=(awaken-sleep).total_seconds(),
    request_to_hibernation_wall_seconds=(sleep-start).total_seconds(),
    deadline_occurred_during_hibernation=sleep<deadline<awaken,
    after_resume='budget guards ended unfinished calculation and waiting queues; no new heavy calculation started',
    no_global_power_policy_changed=True,no_new_automation_scheduled=True,old_si7_automation_remains_paused=True,
    known_calculation_and_waiter_processes_absent_verified=True,
    record_written_UTC=datetime.now(timezone.utc).isoformat(),events=events)
dump(stage/'session_interruption.json',interruption)
with np.load(root/'dense_loading8/partial.npz') as raw:
    partial=raw['hessian_rows'];basis=raw['basis'];gradient=raw['gradient']
    leading=partial[:,:len(partial)];asym=float(np.max(abs(leading-leading.T)))
    values=np.linalg.eigvalsh((leading+leading.T)/2)
    audit=dict(complete_partial_file_integrity_audit=True,rows=partial.shape[0],full_dimension=partial.shape[1],
        raw_all_finite=bool(np.all(np.isfinite(partial))),basis_shape=list(basis.shape),
        known_leading_block_asymmetry_max_eV_A2=asym,
        known_leading_principal_block_minimum_curvature_eV_A2=float(values[0]),
        reduced_gradient_max_eV_A=float(np.max(abs(gradient))),
        full_matrix_completed=False,full_matrix_stability_certified=False,
        interpretation='positive leading principal block does not certify the missing full matrix; no independent derivative checks at this state completed',
        analysis_after_resume=True,new_potential_calls=0,new_DFT=0,new_MD=0,
        partial_file_sha256=hashlib.sha256((root/'dense_loading8/partial.npz').read_bytes()).hexdigest())
dump(stage/'partial_hessian_audit.json',audit)
entries=[]
for label,folder in [('four_row_probe','hessian_probe'),('force100_full','dense_force100'),('loading8_partial','dense_loading8')]:
    value=json.loads((root/folder/'summary.json').read_text(encoding='utf-8'))
    entries.append(dict(name=label,result=folder+'/summary.json',**value))
entries.extend([
    dict(name='probe_standard_calculator_difference',status='complete',new_model_calls=4,result='hessian_probe/standard_calculator_replay.json'),
    dict(name='radial_delta_fit',status='complete_not_adopted',actual_linear_fits=6,new_model_calls=0),
    dict(name='radial_angular_delta_fit',status='complete_not_adopted',actual_linear_fits=6,new_model_calls=0),
    dict(name='structural_descriptors',status='complete',stored_states_reanalyzed=7,new_model_calls=0),
    dict(name='charge_validator',status='reproduced_fixed_and_tested',old_max_probability_loss=0.06866588757887537,new_model_calls=0),
    dict(name='v9_charge_regression',status='actually_rerun',byte_equal_reference_files=5,new_model_calls=0),
    dict(name='hard_charge_boundary',status='complete_synthetic',matrices=24,time_points=120,exact_boundary_controls=27,new_model_calls=0),
    dict(name='charge_spatial_continuum',status='complete_synthetic',sparse_backward_systems=120,analytic_boundary_checks=12,full_vector_replay_states=5190,new_model_calls=0),
    dict(name='implementation_tests',status='complete',tests_passed=63,warnings=2,elapsed_seconds=5.16,source_files_bound=14),
    dict(name='initial_windows_powershell_launcher',status='failed_then_recovered',reason='execution policy; no atomic calculation started in failed attempt'),
    dict(name='first_test_record_wrapper',status='recording_failed_then_recovered',reason='warning-bearing pytest summary parser; tests themselves passed63; final rerun binds source hashes')])
for name in ('dense_return8','dense_loading10','grip_augmented','large_oxide_321_1200',
             'large_oxide_replay','large_delta_transfer','mpa_stratified_comparison','mpa_independent_replay',
             'return_zero_force','zero_force_return_replay','hessian_localization','local_harmonic_diagnostic',
             'local_anharmonic_audit','mpa_prism_transfer'):
    if (root/name).exists():raise ValueError('unexpected final output requires active review: '+name)
    entries.append(dict(name=name,status='not_executed',new_model_calls=0))
ledger=dict(status='final accounting after interrupted eight-hour plan',requested_eight_hours_completed=False,
    actual_new_potential_forward_evaluations=55,actual_autograd_Hessian_rows=1192,
    actual_new_DFT=0,actual_new_MD=0,actual_linear_fits=12,
    full_constrained_Hessians=1,planned_full_constrained_Hessians=4,
    partial_loading8_rows=539,partial_loading8_dimension=648,
    small_probe_rows=4,small_probe_dimension=649,
    potential_call_accounting='1 probe forward +4 independent standard calculator calls +49 completed force100 forwards +1 partial loading8 forward',
    Hessian_row_accounting='4 probe +649 completed force100 +539 partial loading8',
    elapsed_is_not_CPU_time='loading8 elapsed includes Windows hibernation; never summed as continuous compute time',
    large_structure_priority_revision='prepared before prediction: prioritize69 oxygen-containing larger frames;99 larger pure-Si frames deferred; neither executed',
    source_preparation_is_not_execution=True,entries=entries)
dump(stage/'execution_ledger.json',ledger)
print(json.dumps(dict(hibernation_seconds=interruption['reported_hibernation_wall_seconds'],partial=audit,ledger_counts={k:v for k,v in ledger.items() if k!='entries'}),indent=2))
