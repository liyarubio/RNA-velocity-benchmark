"""Verify the curated 1019 code/data release using only the Python standard library."""
from pathlib import Path
from urllib.parse import unquote,urlsplit
import csv,hashlib,json,re

ROOT=Path(__file__).resolve().parents[1];RELEASE=ROOT/'Manuscript_1019'
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()
def main():
    manifest=json.loads((RELEASE/'release_manifest.json').read_text());errors=[]
    files=manifest['copied_source_files']
    for rel,digest in files.items():
        p=RELEASE/rel
        if not p.is_file() or sha(p)!=digest:errors.append(['source_identity',rel])
    index=list(csv.DictReader((RELEASE/'Figure_materials/figure_index.csv').open()))
    expected=[f'Figure{n:02d}' for n in range(1,8)]+[f'ED{n:02d}' for n in range(1,13)]
    if [r['figure'] for r in index]!=expected:errors.append(['figure_order'])
    inventory=list(csv.DictReader((RELEASE/'Figure_materials/file_inventory.csv').open()))
    for row in inventory:
        if row['source_path'] not in files or files[row['source_path']]!=row['sha256']:errors.append(['inventory_binding',row['source_path']])
    links=0
    for p in [ROOT/'README.md']+list(RELEASE.rglob('*.md')):
        # Historical upstream README files retain their original external-workspace references.
        if 'upstream' in p.relative_to(ROOT).parts:continue
        for value in re.findall(r'!?\[[^\]]*\]\((<[^>]+>|[^)\s]+)\)',p.read_text()):
            value=value.strip().strip('<>');u=urlsplit(value)
            if u.scheme or u.netloc or not u.path:continue
            target=(p.parent/unquote(u.path)).resolve();links+=1
            if target==RELEASE/'release_validation.json':continue
            if not target.exists() or not target.is_relative_to(ROOT):errors.append(['local_link',str(p.relative_to(ROOT)),value])
    protected=[]
    for p in RELEASE.rglob('*'):
        if not p.is_file():continue
        if p.stat().st_size>=100_000_000:errors.append(['oversize',str(p.relative_to(ROOT))])
        if p.name in ['.env','id_ed25519','id_rsa'] or p.suffix in ['.pem','.key']:errors.append(['unexpected_credential_file',str(p.relative_to(ROOT))])
        if p.suffix in ['.py','.json','.md','.txt','.csv']:
            text=p.read_text(errors='replace')
            if re.search(r'-----BEGIN (?:OPENSSH |RSA |EC )?PRIVATE KEY-----|gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{30,}',text):errors.append(['possible_secret',str(p.relative_to(ROOT))])
    result={'status':'pass' if not errors else 'fail','figures':len(index),'source_files_checked':len(files),
        'figure_file_mappings':len(inventory),'local_markdown_links_checked':links,'errors':errors,
        'models_rerun':False,'manuscript_included':False,'scope':'Copied-source SHA256, figure order, links, file size and credential-pattern checks.'}
    (RELEASE/'release_validation.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False));return bool(errors)
if __name__=='__main__':raise SystemExit(main())
