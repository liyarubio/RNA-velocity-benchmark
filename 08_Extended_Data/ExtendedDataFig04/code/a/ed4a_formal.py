"""Render ED04a with Mean only and aligned top-three ranks."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / "12_Provenance/current_sources.json").is_file())
sys.path.insert(0, str(ROOT / "09_Core/visualization/mean_top3"))
from mean_top3 import main

if __name__ == "__main__":
    main("ED04a")
