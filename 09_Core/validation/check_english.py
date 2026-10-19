"""Check release text, decoded notebook strings, workbook XML and PDF text for Han characters."""
from pathlib import Path
import argparse
import json
import re
import xml.etree.ElementTree as ET
import zipfile

ROOT = Path(__file__).resolve().parents[2]
SKIP = {'.git', '.venv', '__pycache__', 'demo_output', 'panel_checks'}
HAN = re.compile('[\u3400-\u9fff\uf900-\ufaff\U00020000-\U000323af]')
TEXT = {'.md', '.txt', '.csv', '.json', '.ipynb', '.py', '.R', '.r', '.sh', '.html', '.svg', '.yml', '.yaml', '.toml', '.atr', '.cdt', '.gtr'}


def strings(value):
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for key, child in value.items():
            yield key
            yield from strings(child)
    elif isinstance(value, list):
        for child in value:
            yield from strings(child)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pdf', action='store_true', help='Also extract PDF text; requires PyMuPDF.')
    parser.add_argument('--out', type=Path, help='Optional JSON report path.')
    args = parser.parse_args()
    if args.pdf:
        import fitz
    findings = []
    counts = dict(paths=0, text_files=0, decoded_json_files=0, workbooks=0, pdf_files=0)

    def check(value, location):
        if HAN.search(value):
            findings.append(location)

    for path in sorted(ROOT.rglob('*')):
        rel = path.relative_to(ROOT)
        if SKIP.intersection(rel.parts) or not path.is_file():
            continue
        counts['paths'] += 1
        check(rel.as_posix(), rel.as_posix() + ': filename')
        if path.suffix in TEXT or path.name in {'VERSION', '.gitignore'}:
            content = path.read_text()
            counts['text_files'] += 1
            check(content, rel.as_posix() + ': text')
            if path.suffix in {'.json', '.ipynb'}:
                counts['decoded_json_files'] += 1
                for value in strings(json.loads(content)):
                    check(value, rel.as_posix() + ': decoded JSON')
        elif path.suffix == '.xlsx':
            counts['workbooks'] += 1
            with zipfile.ZipFile(path) as archive:
                for name in archive.namelist():
                    if name.endswith(('.xml', '.rels')):
                        for elem in ET.fromstring(archive.read(name)).iter():
                            for value in [elem.text or '', elem.tail or '', *elem.attrib.values()]:
                                check(value, rel.as_posix() + ': ' + name)
        elif path.suffix == '.pdf' and args.pdf:
            counts['pdf_files'] += 1
            with fitz.open(path) as document:
                for i, page in enumerate(document):
                    check(page.get_text(), f'{rel}: PDF page {i + 1}')
    report = dict(status='pass' if not findings else 'fail', counts=counts,
                  findings=sorted(set(findings)), pdf_text_checked=args.pdf,
                  scope='Current working tree; excludes Git history and local generated outputs. Raster images and vector outlines are not OCR-scanned.')
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
    return bool(findings)


if __name__ == '__main__':
    raise SystemExit(main())
