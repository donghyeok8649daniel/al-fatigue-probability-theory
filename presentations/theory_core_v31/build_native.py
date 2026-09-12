"""Build editable theory slides with installed PowerPoint's document API.

Fallback for a host without @oai/artifact-tool. No GUI automation, execution
policy changes, third-party slide generator, or modification of a source deck.
Outputs go to a fresh private build directory; publication is a separate step.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import io
import zipfile
import xml.etree.ElementTree as ET


POWERSHELL = r'''
$ErrorActionPreference='Stop'
$taskData=Get-Content -LiteralPath $env:THEORY_SOURCE -Raw -Encoding utf8 | ConvertFrom-Json
$taskOut=$env:THEORY_OUTPUT
if(Test-Path -LiteralPath (Join-Path $taskOut 'theory_core_kr_v31.pptx')) {throw 'Refusing to overwrite a built deck'}
[void](New-Item -ItemType Directory -Path $taskOut -Force)
$taskApp=New-Object -ComObject PowerPoint.Application
$taskPres=$taskApp.Presentations.Add(0)
$taskPres.PageSetup.SlideWidth=960
$taskPres.PageSetup.SlideHeight=540
$taskChecks=[System.Collections.Generic.List[object]]::new()
function AddText($slide,$text,$x,$y,$w,$h,$size,$math,$color) {
    $box=$slide.Shapes.AddTextbox(1,$x,$y,$w,$h)
    $box.TextFrame.MarginLeft=0; $box.TextFrame.MarginRight=0
    $box.TextFrame.MarginTop=0; $box.TextFrame.MarginBottom=0
    $box.TextFrame.WordWrap=0
    $box.TextFrame.AutoSize=0
    $matches=[regex]::Matches($text,'([_^])\{([^{}]+)\}')
    $plain=[regex]::Replace($text,'([_^])\{([^{}]+)\}','$2')
    $range=$box.TextFrame.TextRange
    $range.Text=$plain
    $range.Font.Name='Malgun Gothic'
    $range.Font.NameFarEast='Malgun Gothic'
    if($math){$range.Font.Name='Cambria Math'}
    $range.Font.Size=$size
    $range.Font.Color.RGB=$color
    $range.ParagraphFormat.SpaceBefore=0; $range.ParagraphFormat.SpaceAfter=0
    $removed=0
    foreach($match in $matches){
        $start=$match.Index-$removed+1
        $part=$range.Characters($start,$match.Groups[2].Value.Length)
        $part.Font.Size=[single]($size*.70)
        if($match.Groups[1].Value -eq '_'){$part.Font.BaselineOffset=[single](-.22)}
        else{$part.Font.BaselineOffset=[single](.35)}
        $removed+=3
    }
    $taskChecks.Add([pscustomobject]@{slide=$slide.SlideIndex;text=$plain;width=$w;height=$h;bound_width=$range.BoundWidth;bound_height=$range.BoundHeight;font_size=$size;overflow=($range.BoundWidth -gt $w+1 -or $range.BoundHeight -gt $h+1)})
    return $box
}
try {
    $index=0
    foreach($item in $taskData.slides){
        $index++
        $slide=$taskPres.Slides.Add($index,12)
        $slide.FollowMasterBackground=0
        $slide.Background.Fill.Solid()
        $slide.Background.Fill.ForeColor.RGB=16777215
        if($item.cover){
            [void](AddText $slide $item.title 56 135 850 135 42 $false 3154195)
            $y=310
            foreach($line in $item.body){[void](AddText $slide $line 58 $y 842 45 23 $false 6316128);$y+=47}
        } else {
            [void](AddText $slide $item.title 54 42 858 60 32 $false 3154195)
            $eqs=@($item.equations)
            if($null -eq $item.equations){$eqs=@()}
            $y=146
            foreach($line in $eqs){
                $clean=[regex]::Replace($line,'([_^])\{([^{}]+)\}','$2')
                $size=25
                if($clean.Length -gt 87){$size=22}
                elseif($clean.Length -gt 75){$size=23}
                [void](AddText $slide $line 58 $y 844 53 $size $true 7368960)
                $y+=68
            }
            if($eqs.Count){$y=[math]::Max($y+28,372)}else{$y=151}
            foreach($line in $item.body){
                [void](AddText $slide $line 58 $y 848 43 20 $false 4605510)
                $y+=57
            }
        }
        [void](AddText $slide ('LJ / BESSEL · THEORY     '+$index.ToString('00')+' / '+$taskData.slides.Count) 58 506 844 20 11 $false 8421504)
        $notes=$item.notes+"`r`n`r`n출처`r`n"+($item.sources -join "`r`n")
        $slide.NotesPage.Shapes.Placeholders.Item(2).TextFrame.TextRange.Text=$notes
    }
    $taskChecks | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $taskOut 'layout_checks.json') -Encoding utf8
    $taskPres.SaveAs((Join-Path $taskOut 'theory_core_kr_v31.pptx'),24)
    foreach($slide in $taskPres.Slides){$slide.Export((Join-Path $taskOut ('slide_'+$slide.SlideIndex.ToString('00')+'.png')),'PNG',1600,900)}
    $taskPres.SaveAs((Join-Path $taskOut 'theory_core_kr_v31.pdf'),32)
    Write-Output ('Built '+$taskPres.Slides.Count+' slides')
} finally {
    $taskPres.Close()
    if($taskApp.Presentations.Count -eq 0){$taskApp.Quit()}
}
'''


def main() -> None:
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--pdf-tools', type=Path, help='optional private PyMuPDF dependency directory')
    args = parser.parse_args()
    source = Path(__file__).with_name('slides.json').resolve()
    data = json.loads(source.read_text(encoding='utf8'))
    root = source.parents[2]
    missing = sorted({s for slide in data['slides'] for s in slide['sources']
                      if not s.startswith('https://') and not (root / s).is_file()})
    if missing:
        raise ValueError(f'Missing source files: {missing}')
    env = dict(os.environ, THEORY_SOURCE=str(source), THEORY_OUTPUT=str(args.out.resolve()))
    subprocess.run(['powershell', '-NoProfile', '-Command', POWERSHELL], env=env, check=True)
    # PowerPoint writes the host Office account into lastModifiedBy even when
    # Author is explicitly set. Keep machine/account metadata out of Git.
    deck = args.out / 'theory_core_kr_v31.pptx'
    data_bytes = deck.read_bytes()
    buffer = io.BytesIO()
    with zipfile.ZipFile(io.BytesIO(data_bytes)) as zin, zipfile.ZipFile(buffer, 'w', zipfile.ZIP_DEFLATED) as zout:
        for entry in zin.infolist():
            contents = zin.read(entry.filename)
            if entry.filename == 'docProps/core.xml':
                xml = ET.fromstring(contents)
                for node in xml:
                    if node.tag.rsplit('}', 1)[-1] in ('creator', 'lastModifiedBy'):
                        node.text = 'Al fatigue research'
                contents = ET.tostring(xml, encoding='utf-8', xml_declaration=True)
            zout.writestr(entry, contents)
    deck.write_bytes(buffer.getvalue())
    import sys
    if args.pdf_tools:
        sys.path.insert(0, str(args.pdf_tools.resolve()))
    import pymupdf
    pdf_path = deck.with_suffix('.pdf')
    with pymupdf.open(pdf_path) as pdf:
        metadata = pdf.metadata
        metadata.update(author='Al fatigue research', title=data['title'])
        pdf.set_metadata(metadata)
        pdf_bytes = pdf.tobytes(garbage=4, deflate=True)
    pdf_path.write_bytes(pdf_bytes)
    notes = [f"# {data['title']}", '', data['subtitle'], '']
    for i, slide in enumerate(data['slides'], 1):
        notes += [f"## {i:02d}. {slide['title'].replace(chr(10), ' ')}", '']
        for eq in slide.get('equations', []):
            notes += [f'    {eq}', '']
        notes += [*slide['body'], '', slide['notes'], '', '출처:', '']
        notes += [f'- [{s}]({s if s.startswith("https://") else "../../"+s})' for s in slide['sources']]
        notes += ['']
    (args.out / 'theory_notes.md').write_text('\n'.join(notes), encoding='utf8')


if __name__ == '__main__':
    main()
