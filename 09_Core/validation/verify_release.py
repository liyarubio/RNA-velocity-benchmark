"""Verify migrated source identity, figure order and current navigation (stdlib only)."""
from pathlib import Path
from urllib.parse import unquote, urlsplit
import csv
import hashlib
import json
import re

ROOT = Path(__file__).resolve().parents[2]
PROVENANCE = ROOT / '12_Provenance'
LOCAL_OUTPUTS = {'.git', '.venv', '__pycache__', 'demo_output', 'panel_checks'}


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    errors = []
    manifest = json.loads((PROVENANCE / 'path_migration.json').read_text())
    migrations = manifest['files']
    destinations = {}
    for old, entry in migrations.items():
        rel = entry['new_path']
        path = (ROOT / rel).resolve()
        if not path.is_relative_to(ROOT) or not path.is_file():
            errors.append(['missing_migrated_file', old, rel])
            continue
        if rel in destinations:
            errors.append(['duplicate_destination', rel])
        destinations[rel] = old
        if digest(path) != entry['sha256_after']:
            errors.append(['file_identity', rel])
        if entry['change'] == 'moved' and entry['sha256_before'] != entry['sha256_after']:
            errors.append(['pure_move_changed_bytes', old])
        if entry['change'] == 'adapted' and old not in {
            'Manuscript_1019/Code/reproduce.py',
            'Manuscript_1019/Code/plot_saved_panels.py',
            'Manuscript_1019/.gitignore',
        }:
            errors.append(['unexpected_source_edit', old])
    table = list(csv.DictReader((PROVENANCE / 'path_migration.csv').open()))
    if {r['old_path']: {k: v for k, v in r.items() if k != 'old_path'} for r in table} != migrations:
        errors.append(['migration_table_disagrees'])
    index = list(csv.DictReader((PROVENANCE / 'figure_index.csv').open()))
    expected = [f'Figure{i:02}' for i in range(1, 8)] + [f'ED{i:02}' for i in range(1, 13)]
    if [r['figure'] for r in index] != expected:
        errors.append(['figure_order'])
    for row in index:
        if not (ROOT / row['directory'] / 'README.md').is_file():
            errors.append(['figure_entrypoint', row['figure']])
        if not any((ROOT / row['directory'] / 'figures').glob('*.pdf')):
            errors.append(['missing_formal_figure', row['figure']])
    inventory = list(csv.DictReader((PROVENANCE / 'file_inventory.csv').open()))
    for row in inventory:
        old = 'Manuscript_1019/' + row['original_source_path']
        if old not in migrations or migrations[old]['new_path'] != row['source_path']:
            errors.append(['inventory_binding', row['source_path']])
        elif migrations[old]['sha256_after'] != row['sha256']:
            errors.append(['inventory_hash', row['source_path']])
    for row in csv.DictReader((PROVENANCE / 'workbook_sheets.csv').open()):
        old = 'Manuscript_1019/' + row['original_source_path']
        if old not in migrations or migrations[old]['new_path'] != row['source_path']:
            errors.append(['workbook_binding', row['source_path']])
    inherited_prefixes = ('13_Archive/', '12_Provenance/reports/', '12_Provenance/baseline_1019/')
    links = 0
    for path in ROOT.rglob('*.md'):
        rel = path.relative_to(ROOT).as_posix()
        if LOCAL_OUTPUTS.intersection(path.relative_to(ROOT).parts) or rel.startswith(inherited_prefixes):
            continue
        for value in re.findall(r'!?\[[^\]]*\]\((<[^>]+>|[^)\s]+)\)', path.read_text()):
            uri = urlsplit(value.strip('<>'))
            if uri.scheme or uri.netloc or not uri.path:
                continue
            target = (path.parent / unquote(uri.path)).resolve()
            links += 1
            if target == PROVENANCE / 'restructure_validation.json':
                continue
            if not target.is_relative_to(ROOT) or not target.exists():
                errors.append(['local_link', rel, value])
    expected_roots = {f'{i:02}_Figure{i:02}' for i in range(1, 8)} | {
        '08_Extended_Data', '09_Core', '10_Shared', '11_Tables', '12_Provenance', '13_Archive'}
    observed = {p.name for p in ROOT.iterdir() if p.is_dir() and not p.name.startswith('.') and p.name not in LOCAL_OUTPUTS}
    if observed != expected_roots:
        errors.append(['root_directories', sorted(observed ^ expected_roots)])
    secret_pattern = re.compile(r'-----BEGIN (?:OPENSSH |RSA |EC )?PRIVATE KEY-----|gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{30,}')
    for path in ROOT.rglob('*'):
        if LOCAL_OUTPUTS.intersection(path.relative_to(ROOT).parts) or not path.is_file():
            continue
        rel = path.relative_to(ROOT).as_posix()
        if path.is_symlink():
            errors.append(['unexpected_symlink', rel])
        if path.stat().st_size >= 100_000_000:
            errors.append(['oversize', rel])
        if path.suffix in {'.py', '.json', '.md', '.txt', '.csv', '.ipynb', '.sh', '.R'}:
            if secret_pattern.search(path.read_text(errors='replace')):
                errors.append(['possible_secret', rel])
    report = dict(status='pass' if not errors else 'fail', baseline_commit=manifest['baseline_commit'],
                  migrated_files=len(migrations), unchanged_files=sum(r['change'] == 'moved' for r in migrations.values()),
                  figures=len(index), figure_file_mappings=len(inventory), current_local_links_checked=links,
                  inherited_navigation_excluded=list(inherited_prefixes), errors=errors, models_rerun=False)
    (PROVENANCE / 'restructure_validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
    return bool(errors)


if __name__ == '__main__':
    raise SystemExit(main())
