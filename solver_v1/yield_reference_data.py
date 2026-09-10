"""Primary-source tensile YS records; no silent proof-strain or texture mapping."""
import hashlib
import re
import xml.etree.ElementTree as ET


PIGATO_URL = 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13027578/fullTextXML'
PIGATO_DOI = '10.3390/ma19061195'


def parse_pigato_yield_xml(raw):
    """Extract the twelve reported YS numbers, preserving conditions/unknowns.

    The article says YS and cites ISO6892; the retrieved text does not explicitly
    specify a numerical offset. It is NOT relabeled as Rp0.2. Reported uncertainty
    is not reconstructed for specimens with only one test. The three starting
    microstructures are separate records, not one pure-Al mean strength.
    """
    root = ET.fromstring(raw)
    text = lambda node: ''.join(node.itertext()).strip() if node is not None else ''
    dois = [text(n) for n in root.findall('.//article-id') if n.get('pub-id-type') == 'doi']
    if PIGATO_DOI not in dois:
        raise ValueError('unexpected primary-source DOI')
    tables = [n for n in root.findall('.//table-wrap') if text(n.find('label')) == 'Table 2']
    if len(tables) != 1:
        raise ValueError('exactly one source Table 2 required')
    table = tables[0]
    headings = [text(n) for n in table.findall('./table/thead/tr/th')]
    if headings[-3:] != ['6N', '5N5', '5N']:
        raise ValueError('source specimen columns changed')
    rows = table.findall('./table/tbody/tr')
    first = [text(n) for n in rows[0]]
    if re.sub(r'\s+', '', first[0]) != 'YS(MPa)' or rows[0][0].get('rowspan') != '4':
        raise ValueError('source YS units or temperature-block structure changed')
    conditions = {
        '6N': (99.9999, 'fully recrystallized equiaxed; reported Sections3.3/4'),
        '5N5': (99.9995, 'deformation-fragmented substructure; reported Sections3.3/4'),
        '5N': (99.999, 'recovered/subgrain microstructure; reported Section4'),
    }
    records = []
    for index, row in enumerate(rows[:4]):
        cells = [text(n).replace('\u2212', '-') for n in row]
        if index == 0:
            cells = cells[1:]
        if len(cells) != 4:
            raise ValueError('unexpected YS row structure')
        temperature_c = float(cells[0])
        for grade, cell in zip(headings[-3:], cells[1:]):
            match = re.fullmatch(r'(\d+(?:\.\d+)?)\s*(?:±\s*(\d+(?:\.\d+)?))?', cell)
            if match is None:
                raise ValueError('cannot parse source value without inventing a number')
            value = float(match[1]); uncertainty = None if match[2] is None else float(match[2])
            purity, state = conditions[grade]
            records.append(dict(id=f'Pigato2026_{grade}_{temperature_c:g}C_YS',
                material='Al', purity_wt_percent=purity, microstructure=state,
                temperature_C=temperature_c, temperature_K=temperature_c+273.15,
                reported_tensile_YS_MPa=value, reported_uncertainty_MPa=uncertainty,
                uncertainty_convention='source95% confidence; null where only one test reported',
                observable='reported tensile YS; numerical offset not explicitly stated in retrieved text',
                plastic_strain_criterion=None, stress_component='axial',
                loading_direction='parallel to plate rolling direction; no specimen-specific cubic orientation',
                strain_rate_per_second=.00025, strain_rate_relative_tolerance=.20,
                source_pin_spacing_m=None, mobile_dislocation_density_m2=None,
                source_population=None, used_in_fit=False,
                comparable_to_current_source_model=False, source_doi=PIGATO_DOI,
                source_table='Table2', source_method='experiment',
                source_title=text(root.find('.//article-title')), source_year=2026,
                source_authors='Pigato et al.', xml_sha256=hashlib.sha256(raw).hexdigest()))
    return records
