"""Verify the current release hashes, source mapping, figure order and navigation."""
from pathlib import Path
from urllib.parse import unquote,urlsplit
import argparse,csv,hashlib,json,re
ROOT=Path(__file__).resolve().parents[2]
P=ROOT/'12_Provenance'
SKIP={'.git','.venv','__pycache__','demo_output','panel_checks'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path)
    args=parser.parse_args();errors=[]
    hashes=json.loads((P/'current_checksums.json').read_text())
    for rel,h in hashes.items():
        p=(ROOT/rel).resolve()
        if not p.is_relative_to(ROOT) or not p.is_file() or sha(p)!=h:errors.append(['hash',rel])
    mapping=json.loads((P/'current_sources.json').read_text())
    for rel,target in mapping.items():
        p=(ROOT/target).resolve()
        if not p.is_relative_to(ROOT) or not p.is_file():errors.append(['source_path',rel,target])
    source=json.loads((P/'release_1021_sources.json').read_text())
    for rel,e in source['files'].items():
        if mapping.get(rel)!=e['path'] or hashes.get(e['path'])!=e['sha256']:errors.append(['source_binding',rel])
    idx=list(csv.DictReader((P/'figure_index.csv').open()))
    if [x['figure'] for x in idx]!=[f'Figure{i:02}' for i in range(1,8)]+[f'ED{i:02}' for i in range(1,14)]:errors.append(['figure_order'])
    for row in idx:
        if not (ROOT/row['directory']/'README.md').is_file():errors.append(['figure_readme',row['figure']])
        name=row['figure'] if row['figure'].startswith('Figure') else 'Extended_Data_Fig'+row['figure'][2:]
        if not (ROOT/row['directory']/'figures'/(name+'.pdf')).is_file():errors.append(['formal_pdf',row['figure']])
    inventory=list(csv.DictReader((P/'file_inventory.csv').open()))
    for row in inventory:
        if mapping.get(row['original_source_path'])!=row['source_path'] or hashes.get(row['source_path'])!=row['sha256']:errors.append(['figure_binding',row['source_path']])
    for row in csv.DictReader((P/'workbook_sheets.csv').open()):
        if mapping.get(row['original_source_path'])!=row['source_path']:errors.append(['workbook_binding',row['source_path']])
    links=0;secret=re.compile(r'-----BEGIN (?:OPENSSH |RSA |EC )?PRIVATE KEY-----|gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{30,}')
    for path in ROOT.rglob('*'):
        rel=path.relative_to(ROOT).as_posix()
        if SKIP.intersection(path.relative_to(ROOT).parts) or not path.is_file():continue
        if path.is_symlink():errors.append(['symlink',rel])
        if path.stat().st_size>=100_000_000:errors.append(['oversize',rel])
        if path.suffix in {'.py','.md','.json','.txt','.csv','.sh','.ipynb','.R'} and secret.search(path.read_text(errors='replace')):errors.append(['possible_secret',rel])
        if path.suffix!='.md' or rel.startswith(('13_Archive/','12_Provenance/reports/','12_Provenance/baseline_1019/')):continue
        for value in re.findall(r'!?\[[^\]]*\]\((<[^>]+>|[^)\s]+)\)',path.read_text()):
            uri=urlsplit(value.strip('<>'))
            if uri.scheme or uri.netloc or not uri.path:continue
            links+=1;target=(path.parent/unquote(uri.path)).resolve()
            if not target.is_relative_to(ROOT) or not target.exists():errors.append(['link',rel,value])
    report=dict(status='pass' if not errors else 'fail',release='1021-results-public-20261023',checksummed_files=len(hashes),source_bindings=len(source['files']),figure_file_bindings=len(inventory),figures=len(idx),local_links_checked=links,errors=errors,models_rerun=False)
    if args.out:args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2));return bool(errors)
if __name__=='__main__':raise SystemExit(main())
