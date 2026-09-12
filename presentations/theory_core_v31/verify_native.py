"""Structural and text-fit checks plus independent PDF rendering for visual QA."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import sys
import xml.etree.ElementTree as ET
import zipfile


def plain(text):
    return re.sub(r'([_^])\{([^{}]+)\}', r'\2', text)


def compact(text):
    return ''.join(text.split())


def run(build, pdf_tools=None):
    if pdf_tools:
        sys.path.insert(0, str(Path(pdf_tools).resolve()))
    import pymupdf
    import numpy as np
    from PIL import Image

    build = Path(build)
    source = Path(__file__).with_name('slides.json')
    data = json.loads(source.read_text(encoding='utf8'))
    deck = build/'theory_core_kr_v31.pptx'
    pdf_path = deck.with_suffix('.pdf')
    a = '{http://schemas.openxmlformats.org/drawingml/2006/main}'
    with zipfile.ZipFile(deck) as z:
        assert z.testzip() is None
        slides = [n for n in z.namelist() if re.fullmatch(r'ppt/slides/slide\d+.xml', n)]
        assert len(slides) == len(data['slides']) == 40
        for index, slide in enumerate(data['slides'], 1):
            root = ET.fromstring(z.read(f'ppt/slides/slide{index}.xml'))
            text = compact(''.join(e.text or '' for e in root.iter(a+'t')))
            for part in [slide['title'], *slide.get('equations', []), *slide['body']]:
                assert compact(plain(part)) in text, (index, part)
            assert not list(root.iter(a+'blip')), 'equations must stay editable text'
            notes = ET.fromstring(z.read(f'ppt/notesSlides/notesSlide{index}.xml'))
            note_text = ''.join(e.text or '' for e in notes.iter(a+'t'))
            assert slide['notes'] in note_text
            assert all(s in note_text for s in slide['sources'])
        core = ET.fromstring(z.read('docProps/core.xml'))
        assert all(e.text == 'Al fatigue research' for e in core
                   if e.tag.rsplit('}', 1)[-1] in ('creator', 'lastModifiedBy'))
        assert not any(b'C:\\' in z.read(n) for n in z.namelist() if n.endswith('.xml'))
    checks = json.loads((build/'layout_checks.json').read_text(encoding='utf-8-sig'))
    assert not any(r['overflow'] for r in checks)
    images = build/'pdf_render'
    images.mkdir(exist_ok=True)
    rows = []
    with pymupdf.open(pdf_path) as pdf:
        assert len(pdf) == 40
        assert pdf.metadata['author'] == 'Al fatigue research'
        for i, page in enumerate(pdf, 1):
            assert tuple(page.rect)[2:] == (960., 540.)
            assert len(page.get_text()) > 70
            assert '\ufffd' not in page.get_text()
            pix = page.get_pixmap(matrix=pymupdf.Matrix(5/3, 5/3), alpha=False)
            dest = images/f'slide_{i:02}.png'
            pix.save(dest)
            native = np.asarray(Image.open(build/f'slide_{i:02}.png').convert('RGB'), dtype=float)
            independent = np.asarray(Image.open(dest).convert('RGB'), dtype=float)
            rows.append(dict(slide=i, native_pdf_pixel_mean_absolute_difference=float(np.mean(abs(native-independent)))))
    report = dict(slides=40, editable_text=True, notes_and_sources_complete=True,
                  native_text_overflow=False, pdf_independently_rendered=True,
                  # Visual inspection is a separate human/model step, not inferred here.
                  requires_visual_review=True, render_comparison=rows,
                  sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (source, deck, pdf_path)})
    (build/'artifact_checks.json').write_text(json.dumps(report, indent=2), encoding='utf8')
    print(json.dumps({k:v for k,v in report.items() if k != 'render_comparison'}))
    print('maximum PDF/native pixel MAE:', max(r['native_pdf_pixel_mean_absolute_difference'] for r in rows))


if __name__ == '__main__':
    p=argparse.ArgumentParser(__doc__)
    p.add_argument('--build',type=Path,required=True)
    p.add_argument('--pdf-tools',type=Path)
    args=p.parse_args()
    run(args.build,args.pdf_tools)
