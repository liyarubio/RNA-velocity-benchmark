"""Render Figure 3h with methods ordered by descending mean GDC."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / "12_Provenance/current_sources.json").is_file())
sys.path.insert(0, str(ROOT / "09_Core/visualization/figure3h_mean"))
from figure3h_mean import main

if __name__ == "__main__":
    main()
