"""Render the formal ED8a panel using the cell-type colors of ED8c."""
from pathlib import Path
import sys
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'12_Provenance/current_sources.json').is_file())
sys.path.insert(0,str(ROOT/'09_Core/visualization/ed08_celltype'))
from ed08_celltype import main
if __name__=='__main__':main()
