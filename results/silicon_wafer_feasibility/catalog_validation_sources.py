"""Freeze source roles and access limits before new fitting, not fit new data."""
import argparse
import hashlib
import json
from datetime import datetime,timezone,timedelta
from pathlib import Path

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def main(a):
    ROOT=a.root.resolve()
    OUT=a.output
    if OUT.exists(): raise ValueError('fresh source inventory required')
    OUT.mkdir()
    prior=ROOT/'results/silicon_doping_v7/sources/provenance.json'
    info=json.loads(prior.read_text())['datasets']
    records=[]
    for name,value in info.items():
        role={'jaakkola_2014_elastic.csv':'elasticity_temperature_doping',
              'noda_2023_ideal_strength.csv':'static_excess_charge_ideal_strength',
              'chen_2025_wafer.csv':'monotonic_wafer_fracture'}[name]
        records.append(dict(id=name.removesuffix('.csv'),title=value['title'],authors=value['authors'],
            year=value['year'],stable_id=value.get('doi',value.get('stable_id')),
            primary_url=value['url'],role=role,local_numeric_file='results/silicon_doping_v7/sources/'+name,
            local_numeric_sha256=sha(prior.parent/name),previously_examined=True,
            raw_per_specimen_data_obtained=False,first_formation_time_observed=False,
            nonrecovery_test_observed=False,allowed_for_initiation_fit=False,
            existing_provenance_file=str(prior.relative_to(ROOT)).replace('\\','/'),
            access_status=('published tables; authors offer underlying data on request'
                          if role=='monotonic_wafer_fracture' else 'existing published table transcription'),
            warning=value['scope']))
    oxide_path=ROOT/'results/silicon_initiation_v11/thermal_oxide_source.json'
    oxide=json.loads(oxide_path.read_text())
    records.append(dict(id='tsuchiya_2016_oxide',title=oxide['source']['title'],
        authors=oxide['source']['authors'],year=2016,stable_id=oxide['source']['doi'],
        primary_url=oxide['source']['public_repository_url'],role='monotonic_oxide_fracture_stress_sharing',
        local_numeric_file=str(oxide_path.relative_to(ROOT)).replace('\\','/'),local_numeric_sha256=sha(oxide_path),
        previously_examined=True,raw_per_specimen_data_obtained=False,
        first_formation_time_observed=False,nonrecovery_test_observed=False,
        allowed_for_initiation_fit=False,access_status='previously audited user-provided manuscript; no new PDF transcription',
        warning='postmortem origins; inferred flaw sizes not independently measured; oxide thickness convention unresolved'))
    records.append(dict(id='ikehara_2016_orientation',
        title='Crystal orientation-dependent fatigue characteristics in micrometer-sized single-crystal silicon',
        authors=['Tsuyoshi Ikehara','Toshiyuki Tsuchiya'],year=2016,
        stable_id='10.1038/micronano.2016.27',primary_url='https://www.nature.com/articles/micronano201627',
        role='orientation_surface_geometry_final_fatigue_life',previously_examined=True,
        examination_note='publisher full text inspected in this audit; no fitting or numeric digitization performed',
        raw_per_specimen_data_obtained=False,first_formation_time_observed=False,
        nonrecovery_test_observed=False,allowed_for_initiation_fit=False,
        access_status='publisher full text, figures and tables; individual source file not located',
        conditions=dict(total_current_specimens=116,loading='notched resonator bending',
            dimensions_um=[30,10,5],temperature_C=23,relative_humidity_percent=50,
            current_groups=['FR','FT','GU','GR','GS'],previous_groups_excluded=['A','B'],
            ordinate='deflection angle amplitude; FEM conversion needed for stress',
            abscissa='cycles to final failure; CA-equivalent cycles for ramp tests'),
        warning='not direct initiation; growth-law fit is not our initiation law; notch differs from seeded crack; wafer supplies/doping/surface differ; repeated FR curves are not new samples'))
    records.append(dict(id='ando_2016_environment',
        title='Effect of Temperature and Humidity on Degradation of Single-Crystal Silicon Microbeam in MEMS Resonator',
        authors=['Taeko Ando','Mitsuhiro Shikida','Kazuo Sato'],year=2016,
        stable_id='10.18494/SAM.2016.1257',primary_url='https://sensors.myu-group.co.jp/article.php?ss=1257',
        role='temperature_humidity_fatigue_candidate',previously_examined=True,
        examination_note='publisher abstract inspected now; complete PDF fetch failed; numeric figures not transcribed',
        raw_per_specimen_data_obtained=False,first_formation_time_observed=None,
        nonrecovery_test_observed=None,allowed_for_initiation_fit=False,
        access_status='abstract verified; full methods and raw data still needed',
        warning='no numerical thermal/fatigue calibration from abstract; endpoint and stress conversion need full text audit'))
    records.append(dict(id='muhlstein_2001_scs',
        title='High-Cycle Fatigue of Single-Crystal Silicon Thin Films',
        authors=['Christopher L. Muhlstein','Stuart B. Brown','Robert O. Ritchie'],year=2001,
        stable_id='10.1109/84.967383',primary_url='https://www2.lbl.gov/ritchie/Library/PDF/Muhlstein12172001.pdf',
        role='single_crystal_thin_film_total_fatigue_life',previously_examined=True,
        access_status='existing repository bibliography; full-text/raw-figure reaudit not completed this session',
        raw_per_specimen_data_obtained=False,first_formation_time_observed=None,
        nonrecovery_test_observed=None,allowed_for_initiation_fit=False,
        warning='thin film final failure cycles are not first formation; do not transfer environment/size/time to bulk wafer'))
    gap=ROOT/'results/silicon_sw_material_fit/README.md'
    records.append(dict(id='bartok_2018_atomistic',
        title='Machine Learning a General-Purpose Interatomic Potential for Silicon',
        authors=['Albert P. Bartok','James Kermode','Noam Bernstein','Gabor Csanyi'],year=2018,
        stable_id='10.1103/PhysRevX.8.041048',data_id='10.17863/CAM.65004',
        primary_url='https://doi.org/10.17863/CAM.65004',
        primary_article_url='https://doi.org/10.1103/PhysRevX.8.041048',
        role='existing_DFT_energy_force_structure_reference',previously_examined=True,
        access_status='existing atomistic archive and audited predictions; link resolver must be checked before redownload',
        local_provenance_file='results/silicon_sw_material_fit/README.md',local_provenance_sha256=sha(gap),
        raw_per_specimen_data_obtained=False,atomistic_reference_data_obtained=True,
        first_formation_time_observed=False,nonrecovery_test_observed=False,
        allowed_for_initiation_fit=False,
        warning='static/thermal electronic-structure snapshots are not kinetics; previous excluded sets were already seen and are not prospective blind validation'))
    now=datetime.now(timezone.utc)
    ledger=dict(created_utc=now.isoformat(),created_Asia_Seoul=now.astimezone(timezone(timedelta(hours=9))).isoformat(),
        starting_commit='8f9f754c9976b82851a8ad7e2a6a1b365e3bd5c2',sources=records,
        no_new_fit=True,no_new_DFT=True,no_new_MD=True,
        source_builder_sha256=sha(Path(__file__)),
        direct_first_formation_dataset_verified=False,
        direct_nonrecovery_dataset_verified=False,
        matching_nonlinear_equilibrium_or_collective_mobility_dataset_verified=False,
        external_messages_sent=False,
        source_roles_frozen_before_next_fit=True,
        all_sources_inspected_are_no_longer_blind=True,
        required_per_specimen_fields=['specimen_id','parent_wafer_or_batch','wafer_plane','loading_axis',
            'geometry_and_surface_finish','initial_notch_or_crack_state','dopant_species',
            'dopant_chemical_vs_carrier_concentration','oxide_thickness_and_measurement',
            'temperature','humidity_and_oxygen','load_waveform_ratio_frequency','time_or_cycles',
            'first_formation_detection_method_resolution','recovery_load_and_verification_window',
            'endpoint_definition','right_or_interval_censoring','stress_conversion_and_uncertainty'],
        extraction_rule='published summaries, figure digitization and original observations remain distinct; retain unknowns and do not invent specimens',
        fitting_policy='choose new calibration and validation conditions before fitting; never fit failure-life or growth data as first initiation')
    (OUT/'source_ledger.json').write_text(json.dumps(ledger,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(sources=len(records),direct_initiation_data_verified=False,new_fit=False)))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path('.'))
    parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args())
