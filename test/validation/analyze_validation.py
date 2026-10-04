#!/usr/bin/env python3
"""
analyze_validation.py

Utility script for parsing validation CSV outputs, generating summary
comparison markdown and LaTeX tables, and plotting publication-quality figures.

Usage:
    python3 analyze_validation.py --csv <path_to_csv> --out-dir figures/
    python3 analyze_validation.py --latex-table
"""

import os
import sys
import csv
import argparse
from pathlib import Path
from collections import defaultdict

def parse_qr_csv(filepath: Path):
    rows = []
    with open(filepath, "r") as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append({
                "n_sparse_3b": int(r.get("n_sparse_3b", 0)),
                "seed": int(r.get("seed", 0)),
                "variant": r.get("variant", ""),
                "m_total": int(r.get("m_total", 0)),
                "fit_time_s": float(r.get("fit_time_s", 0.0)),
                "nrmse_pct": float(r.get("nrmse_pct", 0.0)),
                "max_rel_pct": float(r.get("max_rel_pct", 0.0)),
                "cosine_sim": float(r.get("cosine_sim", 1.0))
            })
    return rows

def generate_markdown_summary(rows):
    """Summarizes error metrics across variants."""
    grouped = defaultdict(list)
    for r in rows:
        key = (r["variant"], r["n_sparse_3b"])
        grouped[key].append(r)

    lines = []
    lines.append("| Variant | $N_{3b}$ | NRMSE (%) Mean ± Std | Max Rel (%) Mean ± Std | Cosine Sim Mean |")
    lines.append("| :--- | :--- | :--- | :--- | :--- |")

    for (var, n3b), recs in sorted(grouped.items(), key=lambda x: (x[0][0], x[0][1])):
        nrmses = [r["nrmse_pct"] for r in recs]
        max_rels = [r["max_rel_pct"] for r in recs]
        c_sims = [r["cosine_sim"] for r in recs]

        mean_nrmse = sum(nrmses) / len(nrmses)
        std_nrmse = (sum((x - mean_nrmse) ** 2 for x in nrmses) / len(nrmses)) ** 0.5 if len(nrmses) > 1 else 0.0

        mean_rel = sum(max_rels) / len(max_rels)
        std_rel = (sum((x - mean_rel) ** 2 for x in max_rels) / len(max_rels)) ** 0.5 if len(max_rels) > 1 else 0.0

        mean_csim = sum(c_sims) / len(c_sims)

        lines.append(f"| {var} | {n3b} | {mean_nrmse:.4e} ± {std_nrmse:.2e} | {mean_rel:.4e} ± {std_rel:.2e} | {mean_csim:.8f} |")

    return "\n".join(lines)

def generate_latex_table(rows):
    """Generates LaTeX tabular snippet."""
    grouped = defaultdict(list)
    for r in rows:
        key = (r["variant"], r["n_sparse_3b"])
        grouped[key].append(r)

    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"  \centering")
    lines.append(r"  \caption{QR Solver Variants Validation Summary}")
    lines.append(r"  \label{tab:qr_validation}")
    lines.append(r"  \begin{tabular}{llccc}")
    lines.append(r"    \hline\hline")
    lines.append(r"    Variant & $N_{3\mathrm{b}}$ & NRMSE (\%) & Max Rel (\%) & Cosine Sim \\")
    lines.append(r"    \hline")

    for (var, n3b), recs in sorted(grouped.items(), key=lambda x: (x[0][0], x[0][1])):
        nrmses = [r["nrmse_pct"] for r in recs]
        max_rels = [r["max_rel_pct"] for r in recs]
        c_sims = [r["cosine_sim"] for r in recs]
        mean_nrmse = sum(nrmses) / len(nrmses)
        mean_rel = sum(max_rels) / len(max_rels)
        mean_csim = sum(c_sims) / len(c_sims)
        var_clean = var.replace("_", r"\_")
        lines.append(f"    {var_clean} & {n3b} & {mean_nrmse:.2e} & {mean_rel:.2e} & {mean_csim:.7f} \\\\")

    lines.append(r"    \hline\hline")
    lines.append(r"  \end{tabular}")
    lines.append(r"\end{table}")
    return "\n".join(lines)

def plot_scaling_graphs(rows, out_dir: Path):
    try:
        import matplotlib
        if "MPLBACKEND" in os.environ:
            matplotlib.use(os.environ["MPLBACKEND"])
        else:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("[WARN] matplotlib not available; skipping figure generation.")
        return

    out_dir.mkdir(parents=True, exist_ok=True)
    variants = sorted(list(set(r["variant"] for r in rows if r["variant"] != "Full_QR")))
    sparse_levels = sorted(list(set(r["n_sparse_3b"] for r in rows)))

    fig, ax = plt.subplots(figsize=(7, 5), dpi=300)
    for var in variants:
        means = []
        for n3b in sparse_levels:
            subset = [r["nrmse_pct"] for r in rows if r["variant"] == var and r["n_sparse_3b"] == n3b]
            means.append(sum(subset) / len(subset) if subset else 0.0)
        ax.plot(sparse_levels, means, marker="o", lw=1.8, label=var.replace("_", " "))

    ax.set_yscale("log")
    ax.set_xlabel(r"3-Body Sparse Points ($M_{3\mathrm{b}}$)", fontsize=11, fontweight="bold")
    ax.set_ylabel("NRMSE (%) [log scale]", fontsize=11, fontweight="bold")
    ax.set_title("QR Variant Deviation vs Reference Full QR", fontsize=12, fontweight="bold")
    ax.grid(True, which="both", linestyle=":", alpha=0.6)
    ax.legend(frameon=True, fontsize=10)
    fig.tight_layout()

    fig_path = out_dir / "qr_variants_nrmse_scaling.png"
    fig.savefig(fig_path, dpi=300)
    plt.close(fig)
    print(f"[INFO] Saved plot: {fig_path}")

def main():
    parser = argparse.ArgumentParser(description="Analyze JGAP validation outputs")
    parser.add_argument("--csv", type=str, default="", help="Path to input validation CSV")
    parser.add_argument("--out-dir", type=str, default="figures", help="Directory for output figures")
    parser.add_argument("--latex", action="store_true", help="Print LaTeX table code")
    args = parser.parse_args()

    if not args.csv:
        print("Usage: python3 analyze_validation.py --csv <path_to_csv> [--latex] [--out-dir figures]")
        return

    csv_path = Path(args.csv)
    if not csv_path.exists():
        print(f"Error: CSV file not found at {csv_path}")
        sys.exit(1)

    rows = parse_qr_csv(csv_path)
    if not rows:
        print(f"Warning: No rows found in {csv_path}")
        return

    print("\n--- Validation Markdown Summary ---")
    print(generate_markdown_summary(rows))

    if args.latex:
        print("\n--- LaTeX Table Output ---")
        print(generate_latex_table(rows))

    out_dir = Path(args.out_dir)
    plot_scaling_graphs(rows, out_dir)

if __name__ == "__main__":
    main()
