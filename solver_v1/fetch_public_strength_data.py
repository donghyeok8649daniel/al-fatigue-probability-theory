"""Fetch provenance-bound public Al wire data, without executing source notebooks.

Europe PMC exposes the article and supplementary bundle. Source files are
kept in an ignored local cache; a compact manifest is suitable for the repo.
No activation area here is a specimen correlation area or a cell mobility.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import urllib.request
import xml.etree.ElementTree as ET
import zipfile


BASE = 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6265499'


def bounded_fetch(url, maximum_bytes):
    with urllib.request.urlopen(url, timeout=45) as response:
        chunks = []; count = 0
        while True:
            chunk = response.read(min(1024**2, maximum_bytes-count+1))
            if not chunk:
                return b''.join(chunks)
            count += len(chunk)
            if count > maximum_bytes:
                raise ValueError('source download exceeded explicit size budget')
            chunks.append(chunk)


def fetch(directory):
    directory = Path(directory)
    if directory.exists():
        raise FileExistsError('fresh source cache directory required')
    xml = bounded_fetch(BASE+'/fullTextXML', 4*1024**2)
    data = bounded_fetch(BASE+'/supplementaryFiles', 128*1024**2)
    archive = zipfile.ZipFile(io.BytesIO(data))
    listing = []
    for info in archive.infolist():
        listing.append(dict(name=info.filename, bytes=info.file_size,
                            compressed_bytes=info.compress_size))
        if info.filename.endswith('.zip'):
            if info.file_size > 128*1024**2:
                raise ValueError('nested source archive too large for declared audit')
            nested = zipfile.ZipFile(io.BytesIO(archive.read(info)))
            for child in nested.infolist():
                listing.append(dict(name=info.filename+'!'+child.filename,
                    bytes=child.file_size, compressed_bytes=child.compress_size))
    root = ET.fromstring(xml)
    tables = []
    for table in root.findall('.//table-wrap'):
        tables.append(dict(label=''.join(table.find('label').itertext()) if table.find('label') is not None else '',
            caption=''.join(table.find('caption').itertext()) if table.find('caption') is not None else '',
            rows=[[''.join(cell.itertext()) for cell in row] for row in table.findall('.//tr')]))
    directory.mkdir(parents=True)
    (directory/'article.xml').write_bytes(xml)
    (directory/'supplementary.zip').write_bytes(data)
    manifest = dict(doi='10.1016/j.dib.2018.11.047', authors=['S. Verheyden','L. Deillon','A. Mortensen'],
        year=2018, source_url=BASE, xml_sha256=hashlib.sha256(xml).hexdigest(),
        supplementary_sha256=hashlib.sha256(data).hexdigest(),
        supplementary_bytes=len(data), archive_members=listing, tables=tables,
        code_executed=False, production_calibration=False)
    (directory/'manifest.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
    return manifest


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args=parser.parse_args()
    result=fetch(args.out)
    print(json.dumps(dict(supplementary_bytes=result['supplementary_bytes'],
        members=len(result['archive_members']),source_sha256=result['supplementary_sha256']),indent=2))
