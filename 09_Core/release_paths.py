"""Resolve historical release input names to their single current locations."""
from pathlib import Path
import json

REPO_ROOT = Path(__file__).resolve().parents[1]
_MOVES = json.loads((REPO_ROOT / '12_Provenance/path_migration.json').read_text())['files']

def source_path(relative):
    entry = _MOVES['Manuscript_1019/' + str(relative)]
    target = (REPO_ROOT / entry['new_path']).resolve()
    if not target.is_relative_to(REPO_ROOT) or not target.is_file():
        raise FileNotFoundError(target)
    return target
