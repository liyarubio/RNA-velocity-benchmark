import re
import subprocess
import sys
from pathlib import Path


TAGS = ["02", "04", "06", "08", "10"]


def to_ratio(tag: str) -> str:
    # "02" -> "0_2", "10" -> "1_0"
    return f"{tag[0]}_{tag[1]}"


def render_config(original_text: str, tag: str) -> str:
    ratio = to_ratio(tag)
    updates = {
        "count_path": f"Input_data/process_{tag}_count.h5ad",
        "input_path": f"Input_data/process_{tag}.h5ad",
        "output_path": f"result/young2_HSPC_Multi_{ratio}_with_cellinfo/",
        "mv_rna_path": f"Input_data/process_{tag}_RNA.h5ad",
        "mv_atac_path": f"Input_data/process_{tag}_ATAC.h5ad",
        "pv_adata_path": f"Input_data/process_{tag}_lineage.h5ad",
        "pv_pkl_path": f"Input_data/process_{tag}_lineage.pkl",
    }

    text = original_text
    for key, value in updates.items():
        pattern = rf"^{key}\s*=\s*['\"].*?['\"]\s*$"
        repl = f"{key} = '{value}'"
        text, n = re.subn(pattern, repl, text, flags=re.MULTILINE)
        if n != 1:
            raise RuntimeError(f"Failed to update {key} in config.py")
    return text


def required_paths(base_dir: Path, tag: str):
    return [
        base_dir / f"Input_data/process_{tag}.h5ad",
        base_dir / f"Input_data/process_{tag}_count.h5ad",
        base_dir / f"Input_data/process_{tag}_RNA.h5ad",
        base_dir / f"Input_data/process_{tag}_ATAC.h5ad",
        base_dir / f"Input_data/process_{tag}_lineage.h5ad",
        base_dir / f"Input_data/process_{tag}_lineage.pkl",
    ]


def main():
    base_dir = Path(__file__).resolve().parent
    config_path = base_dir / "config.py"
    original = config_path.read_text(encoding="utf-8")

    try:
        for tag in TAGS:
            missing = [p for p in required_paths(base_dir, tag) if not p.exists()]
            if missing:
                print(f"\n[Skip {tag}] Missing inputs:")
                for p in missing:
                    print(f"  - {p}")
                continue

            print(f"\n=== Running velogif for process_{tag} ===")
            config_path.write_text(render_config(original, tag), encoding="utf-8")
            subprocess.run([sys.executable, "Run_velogif.py"], cwd=base_dir, check=True)
            print(f"=== Finished process_{tag} ===")
    finally:
        config_path.write_text(original, encoding="utf-8")
        print("\nRestored original config.py")


if __name__ == "__main__":
    main()
