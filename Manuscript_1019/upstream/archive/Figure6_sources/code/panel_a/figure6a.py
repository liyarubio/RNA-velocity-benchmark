"""Build Figure 6a from local CSV snapshots; no project-global inputs at render time."""
import argparse
from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
METHODS = [
    "velocyto", "scVelo(sto)", "scVelo(dyn)", "VeloAE", "dynamo",
    "VeloVAE", "UniTVelo", "DeepVelo(VAE)", "cellDancer", "veloVI",
    "LatentVelo", "DeepVelo(GCN)", "TFvelo", "TopicVelo", "STT", "biVI",
    "DeepKINET", "SDEvelo", "cell2fate",
]
ACCURACY = ["simulation_cosine", "scLT_cosine", "scNT_cosine",
            "HuBMAP_global", "Dataset10_global", "ReDeeM_global",
            "MPSTA_global", "KPTracer_global"]
NOISE = ["noise_0.4", "noise_1.0", "noise_2.0"]
READS_RATIOS = (0.2, 0.4, 0.6, 0.8)
READS = ["reads_cosine_{:.1f}".format(r) for r in READS_RATIOS]
READS_FILE = "Figure6a_reads_velocity_cosine_similarity.csv"
FILES = {
    "simulation": "Figure6a_simulation_cosine.csv",
    "scLT": "Figure6a_scLT_cosine.csv",
    "scNT": "Figure6a_scNT_cosine.csv",
    "HuBMAP": "Figure6a_HuBMAP_global.csv",
    "Dataset10": "Figure6a_Dataset10_global.csv",
    "ReDeeM": "Figure6a_ReDeeM_global.csv",
    "MPSTA": "Figure6a_MPSTA_global.csv",
    "KPTracer": "Figure6a_KPTracer_global.csv",
    "runtime": "Figure6a_runtime_seconds.csv",
    "noise": "Figure6a_noise_cosine.csv",
}
VELOGIF_CELLS = "Figure6a_VeloGIF_cosine_per_cell.csv"
VELOGIF_METHODS = {"scVelo(sto)", "scVelo(dyn)", "SDEvelo", "TopicVelo",
                   "biVI", "DeepKINET", "cell2fate"}


def summarize_cosine_updates(cells):
    """Validate cell-level reference scores and average valid directions only."""
    required = {"dataset", "method", "cell_id", "cosine", "status", "input_run"}
    if not required.issubset(cells.columns):
        raise ValueError("Missing cell-level cosine metadata")
    cells = cells.copy()
    cells["method"] = cells["method"].map(canonical)
    if cells[["dataset", "method", "cell_id", "status", "input_run"]].isna().any().any():
        raise ValueError("Missing or unknown cosine identifiers")
    if cells.duplicated(["dataset", "method", "cell_id"]).any():
        raise ValueError("Duplicate reference cell within dataset/method")
    if set(cells.dataset) != {"scLT", "scNT"}:
        raise ValueError("Expected scLT and scNT cosine updates")
    known_status = {"valid", "missing_prediction", "zero_vector", "nonfinite_velocity", "coordinate_mismatch"}
    if not set(cells.status).issubset(known_status):
        raise ValueError("Unknown evaluation status")
    cells["cosine"] = pd.to_numeric(cells["cosine"], errors="raise")
    valid = cells.status.eq("valid")
    if not cells.loc[valid, "cosine"].between(-1, 1).all():
        raise ValueError("Valid cosine must be finite and within [-1, 1]")
    if cells.loc[~valid, "cosine"].notna().any():
        raise ValueError("Invalid directions must have blank cosine, not zero")
    rows = []
    for dataset, reference_n in (("scLT", 1756), ("scNT", 1947)):
        subset = cells.loc[cells.dataset.eq(dataset)]
        if set(subset.method) != VELOGIF_METHODS:
            raise ValueError("Unexpected updated methods for " + dataset)
        reference_ids = None
        for method, group in subset.groupby("method", sort=False):
            ids = set(group.cell_id)
            if len(ids) != reference_n or (reference_ids is not None and ids != reference_ids):
                raise ValueError("Reference cells must align across methods in " + dataset)
            reference_ids = ids
            if group.input_run.nunique() != 1:
                raise ValueError("Multiple input runs for one method/dataset")
            scores = group.loc[group.status.eq("valid"), "cosine"]
            rows.append({"dataset": dataset, "method": method,
                         "mean_cosine": scores.mean(), "median_cosine": scores.median(),
                         "reference_cells": reference_n, "valid_cells": len(scores),
                         "invalid_cells": reference_n - len(scores),
                         "valid_fraction": len(scores) / reference_n,
                         "input_run": group.input_run.iloc[0]})
    return pd.DataFrame(rows)


def canonical(value):
    key = re.sub(r"^[0-9]+_", "", str(value)).lower()
    key = re.sub(r"[^a-z0-9]", "", key)
    aliases = {re.sub(r"[^a-z0-9]", "", m.lower()): m for m in METHODS}
    aliases.update({
        "scvelo": "scVelo(sto)", "scvelostochastic": "scVelo(sto)",
        "scvelodynamic": "scVelo(dyn)", "scvelodynamical": "scVelo(dyn)",
        "dynamostochastic": "dynamo", "unitevelo": "UniTVelo",
        "deepvelosa": "DeepVelo(VAE)", "deepvelovaebased": "DeepVelo(VAE)",
        "deepvelogb": "DeepVelo(GCN)", "deepvelogcnbased": "DeepVelo(GCN)",
    })
    return aliases.get(key)


def method_rows(frame):
    frame = frame.copy()
    frame.index = [canonical(x) for x in frame.index]
    frame = frame.loc[frame.index.notna()]
    if not frame.index.is_unique:
        raise ValueError("Duplicate methods after alias normalization")
    return frame.apply(pd.to_numeric, errors="raise").reindex(METHODS)


def rank_index(series, denominator):
    return series.rank(method="average", na_option="keep") / denominator


def summarize_reads_stability(frame):
    """Use Figure 5c's four downsampled comparisons, excluding self-comparison."""
    if set(frame.columns) != {"Method", "Ratio_i", "Ratio_j", "CosineSim"}:
        raise ValueError("Expected Figure 5c pairwise reads cosine columns")
    frame = frame.copy()
    original = frame.Method.copy()
    frame["Method"] = original.map(canonical)
    excluded = original.isin(["MultiVelo", "PhyloVelo"])
    if (frame.Method.isna() & ~excluded).any():
        raise ValueError("Unknown reads stability method")
    for key in ("Ratio_i", "Ratio_j", "CosineSim"):
        frame[key] = pd.to_numeric(frame[key], errors="raise")
    allowed = set(READS_RATIOS) | {1.0}
    if not (frame.Ratio_i.isin(allowed) & frame.Ratio_j.isin(allowed)).all():
        raise ValueError("Unexpected reads ratio")
    if not frame.CosineSim.between(-1, 1).all():
        raise ValueError("Reads cosine must be finite and within [-1, 1]")
    frame = frame.loc[~excluded]
    if frame.duplicated(["Method", "Ratio_i", "Ratio_j"]).any():
        raise ValueError("Duplicate method/reads ratio comparison")
    if not frame.groupby("Method").size().eq(25).all():
        raise ValueError("Each observed method requires all 25 reads comparisons")
    if frame.empty:
        raise ValueError("No Figure 6 methods in reads stability data")
    diagonal = frame.loc[frame.Ratio_i.eq(frame.Ratio_j), "CosineSim"]
    if not np.allclose(diagonal, 1, atol=1e-12, rtol=0):
        raise ValueError("Reads self-comparison must equal 1")
    reference = frame.loc[frame.Ratio_j.eq(1)].pivot(
        index="Method", columns="Ratio_i", values="CosineSim")
    result = reference.reindex(index=METHODS, columns=READS_RATIOS)
    result.columns = READS
    return result


def aggregate(tables, exclude_sclt_scnt=False, cosine_updates=None, reads_stability=None):
    """Average included accuracy domains; retain 60:40 rank aggregation."""
    out = pd.DataFrame(index=pd.Index(METHODS, name="Method"))

    def add_mean(key, frame, skipna=True):
        frame = method_rows(frame)
        if np.isinf(frame.values).any():
            raise ValueError("Infinite input in " + key)
        if ((frame < -1 - 1e-9) | (frame > 1 + 1e-9)).any().any():
            raise ValueError("Signed metric outside [-1, 1]: " + key)
        out[key] = frame.mean(axis=1, skipna=skipna)
        out[key + "_n"] = frame.count(axis=1)

    add_mean(ACCURACY[0], tables["simulation"])
    # The historical R rowMeans for these two sources did not use na.rm.
    if not exclude_sclt_scnt:
        add_mean(ACCURACY[1], tables["scLT"].T, skipna=False)
        add_mean(ACCURACY[2], tables["scNT"].T, skipna=False)
    if cosine_updates is not None:
        if exclude_sclt_scnt:
            raise ValueError("Cannot update excluded scLT/scNT domains")
        # Historical tables have positional row IDs, not matched cell barcodes.
        # Replace whole method-level results; never merge their rows by position.
        for row in summarize_cosine_updates(cosine_updates).itertuples():
            key = row.dataset + "_cosine"
            out.loc[row.method, key] = row.mean_cosine
            out.loc[row.method, key + "_n"] = row.valid_cells
    add_mean(ACCURACY[3], tables["HuBMAP"].T)
    add_mean(ACCURACY[4], tables["Dataset10"].drop(columns=["Median", "Mean"], errors="ignore"))
    add_mean(ACCURACY[5], tables["ReDeeM"])
    for domain in ("MPSTA", "KPTracer"):
        if list(tables[domain].columns) != ["global"]:
            raise ValueError(domain + " must contain only the global metric")
        add_mean(domain + "_global", tables[domain])
    accuracy = [key for key in ACCURACY if key in out]
    out["accuracy_domains_n"] = out[accuracy].count(axis=1)
    out["Accuracy"] = ((out[accuracy] + 1) / 2).mean(axis=1)

    seconds = tables["runtime"].copy()
    seconds.index = seconds.index.astype(int)
    times = method_rows(seconds.loc[[1000, 10000, 50000]].T)
    if (times <= 0).any().any() or np.isinf(times.values).any():
        raise ValueError("Observed runtimes must be finite and positive")
    inv = 1 / times
    # 50K is displayed but was never part of the original scalability score.
    out["Scalability"] = inv[[1000, 10000]].mean(axis=1)
    out["scalability_sizes_n"] = times[[1000, 10000]].count(axis=1)
    lower, upper = np.nanmin(inv.values) - 1, np.nanmax(inv.values)
    for size in (1000, 10000, 50000):
        out["runtime_{}_seconds".format(size)] = times[size]
        out["efficiency_{}".format(size)] = (inv[size] - lower) / (upper - lower)

    if reads_stability is not None:
        reads = summarize_reads_stability(reads_stability)
        out[READS] = reads
        out["stability_ratios_n"] = reads.count(axis=1)
        out["Stability"] = reads.mean(axis=1)
    else:
        noise = method_rows(tables["noise"])
        noise.loc["LatentVelo"] = np.nan  # Explicit historical exclusion only.
        for suffix, key in zip((2, 5, 10), NOISE):
            ratios = []
            for model in ("Mono-kinetic", "Transcriptional_boost", "Multi-forward_branching"):
                baseline = noise[model + "_1"].replace(0, np.nan)
                ratios.append(noise[model + "_" + str(suffix)] / baseline)
            ratios = pd.concat(ratios, axis=1)
            out[key] = ratios.mean(axis=1)
            out[key + "_models_n"] = ratios.count(axis=1)
        out["Stability"] = out[NOISE].mean(axis=1)
    for group in ("Accuracy", "Scalability", "Stability"):
        out[group + "_rank_index"] = rank_index(out[group], len(METHODS))
    # Integer weights avoid splitting mathematically tied scores by roundoff.
    out["Overall_Score"] = (3 * out.Accuracy.rank(method="average") +
                            2 * out.Scalability.rank(method="average")) / (5 * len(METHODS))
    for key in ("Overall_Score", "Accuracy", "Scalability", "Stability"):
        out[key + "_rank"] = out[key].rank(ascending=False, method="min", na_option="keep")
    if np.isinf(out.values).any():
        raise ValueError("Non-finite calculation result")
    return out


def rank_colors(series, palette):
    count = series.count()
    positions = (series.rank(method="average") - 1) / max(count - 1, 1)
    cmap = plt.get_cmap(palette + "_r")
    return [cmap(p) if np.isfinite(p) else (0, 0, 0, 0) for p in positions]


def plot(summary, outdir, stem="Figure6a_overall_performance"):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12,
                         "pdf.fonttype": 42, "svg.fonttype": "none"})
    columns = [
        ("Overall_Score", "Overall Score", "Reds", True),
        ("Accuracy_rank_index", "Accuracy rank index", "Oranges", True),
        (ACCURACY[0], "Cosine similarity in simulation datasets (n=42)", "Oranges", False),
        (ACCURACY[1], "Cosine similarity in scLT dataset", "Oranges", False),
        (ACCURACY[2], "Cosine similarity in scNT dataset", "Oranges", False),
        (ACCURACY[3], "GDC in HuBMAP datasets (n=176)", "Oranges", False),
        (ACCURACY[4], "GDC in classic datasets (n=10)", "Oranges", False),
        (ACCURACY[5], "GDC in ReDeeM dataset", "Oranges", False),
        (ACCURACY[6], "GDC in spatial transcriptomics (MPSTA)", "Oranges", False),
        (ACCURACY[7], "GDC in tumor dataset (KPTracer)", "Oranges", False),
        ("Scalability_rank_index", "Scalability rank index", "Greens", True),
        ("efficiency_1000", "Computational efficiency (1K cells)", "Greens", False),
        ("efficiency_10000", "Computational efficiency (10K cells)", "Greens", False),
        ("efficiency_50000", "Computational efficiency (50K cells)", "Greens", False),
        ("Stability_rank_index", "Stability rank index", "Blues", True),
        (NOISE[0], "Percentage of performance (noise level 0.4)", "Blues", False),
        (NOISE[1], "Percentage of performance (noise level 1.0)", "Blues", False),
        (NOISE[2], "Percentage of performance (noise level 2.0)", "Blues", False),
    ]
    columns.extend((key, "Velocity cosine (reads {:.1f}x vs 1.0x)".format(ratio), "Blues", False)
                   for key, ratio in zip(READS, READS_RATIOS))
    columns = [column for column in columns if column[0] in summary]
    fig, ax = plt.subplots(figsize=(21.6 + max(0, len(columns) - 18) * 1.2
                                    if len(columns) >= 18 else 20, 16))
    fig.subplots_adjust(left=.025, right=.99, bottom=.025, top=.99)
    yvals = -np.arange(len(summary)) * 1.2
    positions, x, previous = [], 5.5, None
    for _, _, palette, bar in columns:
        if previous is not None:
            x += .42 if previous != palette else .22
        width = 1.7 if bar else 1
        positions.append((x, width))
        x += width
        previous = palette
    for i, y in enumerate(yvals):
        if i % 2 == 0:
            ax.add_patch(Rectangle((0, y - .55), x + .2, 1.1, facecolor="#d6d6d6", edgecolor="none"))
        ax.text(.2, y, summary.index[i], va="center", fontsize=14)
    for (key, label, palette, bar), (left, width) in zip(columns, positions):
        values = summary[key]
        if key in ACCURACY or key in READS:
            values = (values + 1) / 2
        colors = rank_colors(values, palette)
        ranks = values.rank(ascending=False, method="min")
        ax.text(left + width / 2, .96, label, rotation=30, rotation_mode="anchor",
                ha="left", va="bottom", fontsize=11.5, fontweight="bold" if bar else "normal")
        ax.plot([left + width / 2] * 2, [.65, .82], color="#666666", linewidth=.65)
        for i, (value, y) in enumerate(zip(values, yvals)):
            if pd.isna(value):
                continue
            if bar:
                ax.add_patch(Rectangle((left, y - .55), width * value, 1.1,
                                       facecolor=colors[i], edgecolor="#707070", linewidth=.9))
                if ranks.iloc[i] <= 6:
                    ax.text(left + width * value / 2, y, str(int(ranks.iloc[i])),
                            ha="center", va="center", fontsize=13)
            elif value >= 0:
                # Ratios slightly greater than 1 are retained, not capped.
                ax.add_patch(Circle((left + width / 2, y), .52 * np.sqrt(value),
                                    facecolor=colors[i], edgecolor="#707070", linewidth=.85))
    bottom = yvals[-1] - 2
    ax.text(10, bottom, "Ranking", ha="center", fontsize=12)
    for j, palette in enumerate(("Reds", "Oranges", "Greens", "Blues")):
        cmap = plt.get_cmap(palette)
        for k in range(5):
            ax.add_patch(Rectangle((8.5 + j * .65, bottom - .5 - (k + 1) * .35), .53, .35,
                                   facecolor=cmap(k / 4), edgecolor="#777777", linewidth=.7))
    ax.annotate("", xy=(8.05, bottom - .45), xytext=(8.05, bottom - 2.25),
                arrowprops={"arrowstyle": "->", "color": "#333333"})
    ax.text(14.7, bottom, "Score", ha="center", fontsize=12)
    for center, value in zip(np.linspace(13, 16.4, 6), (0, .04, .16, .36, .64, 1)):
        ax.add_patch(Circle((center, bottom - 1.2), max(.025, .52 * np.sqrt(value)),
                            facecolor="white", edgecolor="#707070", linewidth=.85))
    ax.text(13, bottom - 2.15, "0%", ha="center", fontsize=11)
    ax.text(16.4, bottom - 2.15, "100%", ha="center", fontsize=11)
    note = "Blank = unavailable. Accuracy uses available domains; coverage is listed in the summary CSV."
    if "_VeloGIF_updated" in stem:
        note += "\nscLT/scNT combine historical and updated runs with different preprocessing; interpret ranks with caution."
    if READS[0] in summary:
        note += ("\nStability = mean velocity cosine at 0.2/0.4/0.6/0.8x reads vs 1.0x (Figure 5c); not included in Overall Score."
                 "\nCosine circle area uses (cosine + 1) / 2; signed values are retained in CSV. TopicVelo*: integrated transitions projected into fixed full-depth gene coordinates.")
    ax.text(0, bottom - 3, note,
            fontsize=10, color="#444444", va="top" if READS[0] in summary else "baseline")
    ax.text(0, 1.15, "a", fontsize=26, fontweight="bold")
    ax.set_xlim(-.15, x + 7.8)
    ax.set_ylim(bottom - (5 if READS[0] in summary else 3.5), 7.2)
    ax.set_aspect("equal")
    ax.axis("off")
    for ext in ("svg", "pdf", "png"):
        fig.savefig(outdir / (stem + "." + ext), dpi=300,
                    bbox_inches="tight", pad_inches=.15, facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "reproduced")
    parser.add_argument("--reads-stability", action="store_true",
                        help="Replace noise stability with Figure 5c reads-downsampling cosine")
    versions = parser.add_mutually_exclusive_group()
    versions.add_argument("--exclude-sclt-scnt", action="store_true",
                        help="Exclude scLT/scNT from accuracy and plot; save separately with a suffix")
    versions.add_argument("--velogif-updates", action="store_true",
                          help="Restore scLT/scNT with local VeloGIF cosine updates; save a separate version")
    args = parser.parse_args()
    tables = {k: pd.read_csv(args.data_dir / name, index_col=0,
                             float_precision="round_trip" if k in ("MPSTA", "KPTracer") else None)
              for k, name in FILES.items()
              if not (args.exclude_sclt_scnt and k in ("scLT", "scNT"))
              and not (args.reads_stability and k == "noise")}
    updates = (pd.read_csv(args.data_dir / VELOGIF_CELLS, float_precision="round_trip")
               if args.velogif_updates else None)
    reads = (pd.read_csv(args.data_dir / READS_FILE, float_precision="round_trip")
             if args.reads_stability else None)
    summary = aggregate(tables, exclude_sclt_scnt=args.exclude_sclt_scnt,
                        cosine_updates=updates, reads_stability=reads)
    stem = "Figure6a_overall_performance"
    if args.exclude_sclt_scnt:
        stem += "_without_scLT_scNT"
    elif args.velogif_updates:
        stem += "_VeloGIF_updated"
    if args.reads_stability:
        stem += "_reads_stability"
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output_dir / (stem + "_summary.csv"), na_rep="NA")
    plot(summary, args.output_dir, stem=stem)
    print(summary[["accuracy_domains_n", "Overall_Score", "Overall_Score_rank"]].to_string())
    print("Outputs:", args.output_dir)


if __name__ == "__main__":
    main()
