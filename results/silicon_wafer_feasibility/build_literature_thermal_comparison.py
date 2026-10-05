"""Traceable published-table comparison and cached local harmonic quantum audit.

No fitting, new potential calls, new DFT/MD, initiation counts or clock mapping.
"""
from __future__ import annotations
import argparse
import csv
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.constants import Boltzmann,Planck,electron_volt,atomic_mass
from solver_v1.silicon_thermal_validity import oscillator_statistics,validate_fixed_cartesian

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def dump(path,value):path.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
def readcsv(path):
    with path.open(encoding='utf-8',newline='') as handle:return list(csv.DictReader(handle))

def main(a):
    if a.output.exists():raise ValueError('fresh output required')
    a.output.mkdir(parents=True)
    ledger=json.loads((a.root/'results/silicon_validation_targets/source_ledger.json').read_text())
    # The old ledger is retained as the previous access record; this is an addendum.
    quotes=[
        ('jaakkola_2014_elastic','are affected to within a few percent by increased doping over the tested wafers.',
         'doping discussion, manuscript PDF p7; Table III p8, Eq8','탄성계수 변화의 크기와 농도·온도 범위를 제한한다.'),
        ('noda_2023_ideal_strength','not doping additional atomic elements/impurities, but purely electron doping.',
         'Methods, published PDF p2; Fig7/8, p5–6','전하만 변화시킨 DFT와 실제 도펀트 원자를 구분한다.'),
        ('chen_2025_wafer','no statistically significant differences in the fracture strengths among the six types of silicon wafers.',
         'Section III.C, Tables I and III','이 시험 조건의 비유의성을 도핑 효과 부재의 보편 법칙으로 쓰지 않는다.'),
        ('tsuchiya_2016_oxide','The nominal strength decreased first with the growth of oxide layer to 50 nm',
         'provided manuscript p5, Fig5 discussion; Tables IV/V p14','산화층의 응력 분담과 공정 후 남는 결함을 함께 검사한다.'),
        ('ikehara_2016_orientation','the first principal stress was a good criterion for fatigue fracture.',
         'Discussion; Figure4 and Tables2/3; Methods','저자들이 비교한 최종 피로 파단과 응력 환산의 범위에서만 인용한다.'),
        ('ando_2016_environment','the fatigue strength was markedly decreased under a highly humid condition',
         'publisher abstract only','습도·온도 상호작용을 검토할 근거이며 정량 보정 표로 쓰지 않는다.'),
        ('muhlstein_2001_high_cycle','observed to decrease monotonically before the specimen finally failed at the notch.',
         'Section III.A, journal pp595–596 (PDF pp3–4), Figures3–5','공진 변화는 선행 손상의 간접 관측이며 첫 균열의 직접 시간 기록이 아니다.'),
        ('bartok_2018_potential','getting the right answer for the wrong reasons.',
         'Section I.A, journal041048-2; Sections III/IV','개별 수치 일치와 에너지·힘·진동·파괴 경로의 공동 검증을 구분한다.')]
    # IDs can be looked up by stable DOI/title if the historical name differs.
    byid={s['id']:s for s in ledger['sources']}
    entries=[]
    for key,quote,locator,role in quotes:
        if key not in byid:
            token={'ando_2016_environment':'1257','muhlstein_2001_high_cycle':'967383','bartok_2018_potential':'041048'}[key]
            matches=[s for s in ledger['sources'] if token in str(s)]
            if len(matches)!=1:raise ValueError('source ID not unique: '+key)
            source=matches[0]
        else:source=byid[key]
        if len(quote.split())>25:raise ValueError('excerpt exceeds 25 words')
        entries.append(dict(id=source['id'],title=source['title'],authors=source['authors'],year=source['year'],
            stable_id=source['stable_id'],primary_url=(
                'https://journals.aps.org/prx/abstract/10.1103/PhysRevX.8.041048'
                if source['stable_id']=='10.1103/PhysRevX.8.041048' else source['primary_url']),excerpt=quote,
            excerpt_word_count=len(quote.split()),locator=locator,research_use=role,
            access_now='abstract only' if key.startswith('ando') else 'relevant full-text passages/figures inspected',
            per_specimen_raw_data_obtained=False,first_formation_time_observed=False,
            nonrecovery_test_observed=False,eligible_for_first_initiation_fit=False))
    dump(a.output/'citation_addendum.json',dict(created_utc=datetime.now(timezone.utc).isoformat(),
        previous_ledger_sha256=sha(a.root/'results/silicon_validation_targets/source_ledger.json'),sources=entries,
        author_contact_performed=False,prior_ledger_modified=False))
    paths=[a.root/'results/silicon_doping_v7/sources'/name for name in
           ('jaakkola_2014_elastic.csv','noda_2023_ideal_strength.csv','chen_2025_wafer.csv')]
    oxide_path=a.root/'results/silicon_initiation_v11/thermal_oxide_source.json'
    j,n,c=[readcsv(p) for p in paths];oxide=json.loads(oxide_path.read_text())
    tables=a.output/'published_tables';tables.mkdir()
    for path in paths+[oxide_path]:(tables/path.name).write_bytes(path.read_bytes())
    fields=['source','sample','observable','value','unit','summary_type','species','carrier_type',
            'carrier_cm3','temperature_C','RH_percent','orientation_or_axis','geometry','loading',
            'sample_count','uncertainty','uncertainty_kind','locator','fit_role','first_initiation_fit_eligible']
    rows=[]
    def add(**kw):
        value={key:'' for key in fields};value.update(kw);value['first_initiation_fit_eligible']='false';rows.append(value)
    for sample in j:
        for key,error in [('c11_GPa',.3),('c12_GPa',.1),('c44_GPa',.2)]:
            add(source='Jaakkola2014',sample=sample['sample_id'],observable=key[:-4].upper(),value=sample[key],unit='GPa',
                summary_type='resonator-inferred cubic stiffness at 25 C',species=sample['species'],
                carrier_type=sample['carrier_type'],carrier_cm3=sample['carrier_nominal_cm3'],temperature_C=25,
                orientation_or_axis='cubic tensor',geometry='MEMS resonators',loading='resonance inverse problem',
                uncertainty=error,uncertainty_kind='reported RMS budget; not independent joint CI',
                locator='Table III; temperature Eq8',fit_role='elasticity validation; carrier inferred from resistivity')
    for sample in n:
        add(source='Noda2023',sample=sample['carrier_type']+'_'+sample['excess_carrier_cm3'],observable='ideal_strength',
            value=sample['ideal_strength_LDA_GPa'],unit='GPa',summary_type='static defect-free DFT maximum stress',
            species='no added dopant atom; compensating background',carrier_type=sample['carrier_type'],
            carrier_cm3=sample['excess_carrier_cm3'],orientation_or_axis='[111]',geometry='periodic perfect crystal',
            loading='homogeneous tensile DFT; LDA',locator='Fig7/8, supplementary S3/S4/S5',
            fit_role='charge response diagnostic; no ordinary-concentration extrapolation')
    for sample in c:
        add(source='Chen2025',sample=sample['sample_id'],observable='Weibull_characteristic_final_strength',
            value=sample['characteristic_strength_GPa'],unit='GPa',summary_type='characteristic strength, not arithmetic mean',
            species=sample['species'],carrier_cm3=sample['carrier_cm3'],orientation_or_axis='(100)',
            geometry='300 mm Cz wafer; 775 um; 17x17 mm slices',loading='ball-on-ring; 0.2 mm/min',
            sample_count=sample['sample_count'],uncertainty=sample['standard_error_GPa'],uncertainty_kind='reported SE',
            locator='Tables I/III; section III.C',fit_role='final monotonic fracture with oxygen ~7.9e17 cm^-3')
    for sample in oxide['rows']:
        common=dict(source='Tsuchiya2016',sample=sample['group'],unit='GPa',temperature_C=26,RH_percent=50,
            orientation_or_axis='(100), <110>',geometry='4x5 um; lengths120/600 um',
            loading='monotonic tensile, stage0.5 um/s',locator='Table IV, supplied PDF p14')
        add(**common,observable='nominal_final_strength',value=sample['nominal_strength_GPa'],
            summary_type='mean nominal stress',fit_role='final fracture and processing-history check')
        for key,label in [('FEM_Si_axial_GPa','Si_local_stress_at_500nm_from_corner'),
                          ('FEM_oxide_axial_GPa','oxide_surface_local_stress')]:
            if sample[key] is not None:
                add(**common,observable=label,value=sample[key],summary_type='FEM output, not phase average or independent measurement',
                    fit_role='same measured load; stress-sharing audit only')
    with (a.output/'numeric_comparison.csv').open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,lineterminator='\n');writer.writeheader();writer.writerows(rows)
    dump(a.output/'table_metadata.json',dict(row_count=len(rows),blank_means_unreported_or_not_applicable=True,
        numeric_inputs_sha256={str(p.relative_to(a.root)).replace('\\','/'):sha(p) for p in paths+[oxide_path]},
        no_new_fit=True,no_new_digitization=True,statistics_must_not_be_interchanged=True,
        missing_numeric_fatigue_raw_sources=['Ikehara2016','Ando2016','Muhlstein2001'],
        missing_initiation_endpoints_for_all_sources=True))

    harmonic=[];mode_arrays={};source_paths=[]
    data=[]
    for name in ('loading8','return8'):
        path=a.root/f'results/silicon_initiation_v13/dense_{name}/raw_hessian.npz';source_paths.append(path)
        with np.load(path) as z:data.append({k:z[k].copy() for k in z.files})
    endpoint=(data[1]['positions']-data[0]['positions']).ravel();endpoint/=np.linalg.norm(endpoint)
    for name,source in zip(('initial8','returned8'),data):
        free=np.any(source['basis'].reshape(len(source['positions']),3,-1)!=0,axis=(1,2))
        h,basis_checks=validate_fixed_cartesian(source['hessian'],source['basis'],free)
        eigenvalues,eigenvectors=np.linalg.eigh(h)
        if np.any(eigenvalues<=0):raise ValueError('positive local modes required; do not clip')
        mass=28.0855*atomic_mass
        frequency_Hz=np.sqrt(eigenvalues*electron_volt/1e-20/mass)/(2*np.pi)
        y=Planck*frequency_Hz/(Boltzmann*300.)
        ratio=(y/2)/np.tanh(y/2)
        var=Boltzmann*300./electron_volt/eigenvalues
        existing=oscillator_statistics(eigenvalues,300.)
        np.testing.assert_allclose(existing['variance_ratio'],ratio,rtol=3e-15,atol=2e-15)
        np.testing.assert_allclose(existing['omega_rad_ps']/(2*np.pi),frequency_Hz/1e12,rtol=3e-15,atol=1e-14)
        projected=source['basis'].T@endpoint
        weights=(projected@eigenvectors)**2
        classical=float(weights@var);quantum=float(weights@(var*ratio))
        direct=Boltzmann*300./electron_volt*float(projected@np.linalg.solve(h,projected))
        np.testing.assert_allclose(direct,classical,rtol=1e-12)
        weights_classical=weights*var/classical
        harmonic.append(dict(name=name,modes=len(eigenvalues),basis_checks=basis_checks,
            minimum_frequency_THz=float(frequency_Hz.min()/1e12),maximum_frequency_THz=float(frequency_Hz.max()/1e12),
            minimum_hnu_over_kBT=float(y.min()),maximum_hnu_over_kBT=float(y.max()),
            maximum_quantum_to_classical_mode_variance=float(ratio.max()),
            endpoint_projection_classical_variance_A2=classical,endpoint_projection_quantum_variance_A2=quantum,
            endpoint_projection_variance_ratio=quantum/classical,
            endpoint_classical_variance_fraction_from_hnu_over_kBT_gt1=float(weights_classical[y>1].sum()),
            total_coordinate_variance_ratio=float((var*ratio).sum()/var.sum()),
            largest_gradient_component_eV_A=float(abs(source['gradient']).max()),
            local_frozen_hessian_approximation=True,endpoint_direction_is_not_verified_crack_coordinate=True,
            anharmonic_quantum_sampling_performed=False,kBT_adjusted=False,physical_clock_calibrated=False))
        for key,value in [('eigenvalues_eV_A2',eigenvalues),('frequency_THz',frequency_Hz/1e12),
                          ('hnu_over_kBT',y),('variance_ratio',ratio),('endpoint_weights',weights),
                          ('endpoint_classical_variance_weights',weights_classical)]:mode_arrays[name+'_'+key]=value
    bulk_path=a.root/'results/silicon_sw_material_fit/bulk_optical_diagnostic.json'
    bulk=json.loads(bulk_path.read_text());nu=bulk['reference_gamma_optical_THz'][0]*1e12
    y=Planck*nu/(Boltzmann*300.)
    dump(a.output/'harmonic_quantum_diagnostic.json',dict(temperature_K=300.,thermal_energy_eV=Boltzmann*300./electron_volt,
        source_sha256={str(p.relative_to(a.root)).replace('\\','/'):sha(p) for p in source_paths+[bulk_path]},
        local_modes=harmonic,bulk_reference=dict(frequency_THz=nu/1e12,hnu_over_kBT=y,
            variance_ratio=(y/2)/np.tanh(y/2),provenance='existing CASTEP DFT benchmark; not experimental measurement',
            applicable_to_current_crack_coordinate=False),
        formula='C_quantum=V diag[(kBT/lambda_j)*(hbar omega_j/(2kBT))*coth(hbar omega_j/(2kBT))] V.T',
        assumptions='stable frozen local harmonic Hessian; full equal-mass orthonormal Cartesian basis; no new eigenvalue clipping',
        primary_theory_url='https://ocw.mit.edu/courses/8-333-statistical-mechanics-i-statistical-mechanics-of-particles-fall-2013/f38da5fe5bb800fb411f1c253ead20a0_b1P0hurY6UE.pdf',
        primary_theory_locator='Mehran Kardar, lecture20, transcript pp1–5; harmonic quantum partition function and mode spectrum',
        new_potential_calls=0,new_MD=0,new_DFT=0,quantum_PMF_estimated=False,material_approved=False,
        energy_quantization_does_not_change_physical_kBT=True,production_changed=False))
    np.savez_compressed(a.output/'harmonic_quantum_modes.npz',**mode_arrays)
    print(json.dumps(dict(numeric_rows=len(rows),harmonic=harmonic),indent=2),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args())
