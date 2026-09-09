"""Optional research-only retrieval of the published experimental benchmark.

PDFs/images stay in ignored .cache; no third-party full text is committed.
Does not add a PDF dependency to the desktop application or numerical solver.
Install pypdf into .cache/research-pdf-tools to use the extraction option.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT/'.cache/strength-reference'
URLS = {
    'article': 'https://eprints.whiterose.ac.uk/id/eprint/126662/1/Article%20Text.pdf',
    # Official ORIGINAL bundle of item d8d8ebd9-c735-4079-a5ce-bbb0cbfffbd1,
    # DOI 10.1038/nmat4911. The legacy /record/228486/files route returns405.
    'supplement': 'https://infoscience.epfl.ch/server/api/core/bitstreams/2cff7b76-19e7-4591-b7ca-370af1bcd090/content',
    'annealed_article': 'https://infoscience.epfl.ch/server/api/core/bitstreams/417397b5-3810-49ce-9e27-7ef484515775/content',
}


def fetch(which):
    CACHE.mkdir(parents=True,exist_ok=True)
    path=CACHE/(which+'.pdf')
    if not path.exists():
        with urlopen(Request(URLS[which],headers={'User-Agent':'research-reference-reader'}),timeout=40) as response:
            data=response.read()
        if not data.startswith(b'%PDF-'):
            raise ValueError('source did not return a PDF; no false cache saved')
        path.write_bytes(data)
    data=path.read_bytes()
    return path,dict(url=URLS[which],sha256=hashlib.sha256(data).hexdigest(),bytes=len(data))


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('which',choices=URLS)
    parser.add_argument('--pages',type=int,nargs='+',default=[5],help='ZERO-based PDF pages')
    parser.add_argument('--images',action='store_true')
    args=parser.parse_args(); path,metadata=fetch(args.which)
    sys.path.insert(0,str(ROOT/'.cache/research-pdf-tools'))
    from pypdf import PdfReader
    reader=PdfReader(path)
    print(json.dumps(metadata),flush=True)
    for number in args.pages:
        page=reader.pages[number]
        print(f'PDF PAGE {number+1}\n'+page.extract_text(),flush=True)
        if args.images:
            for i,img in enumerate(page.images):
                # Lossless original embedded figure, not an AI-generated figure.
                target=CACHE/f'{args.which}_page{number+1}_{i}{Path(img.name).suffix}'
                target.write_bytes(img.data)
                print('embedded_image='+str(target.relative_to(ROOT)),flush=True)
                if target.suffix=='.jp2':
                    # Lossless decoding for viewers without JPEG2000 support;
                    # do not alter, crop, interpolate, or manufacture any pixels.
                    decoded=target.with_suffix('.png')
                    img.image.save(decoded)
                    print('decoded_image='+str(decoded.relative_to(ROOT)),flush=True)


if __name__=='__main__':
    main()
