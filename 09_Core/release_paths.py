"""Resolve original input names to the current, figure-organized release."""
from pathlib import Path
import json
REPO_ROOT = Path(__file__).resolve().parents[1]
_MAPPING = json.loads((REPO_ROOT / '12_Provenance/current_sources.json').read_text())
def source_path(relative):
    target = (REPO_ROOT / _MAPPING[str(relative)]).resolve()
    if not target.is_relative_to(REPO_ROOT) or not target.is_file():
        raise FileNotFoundError(target)
    return target
