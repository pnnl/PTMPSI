from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any
from urllib.parse import unquote
import traceback

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


@dataclass
class ArtifactPaths:
    summary_md: Path
    filtered_tsv: Path
    plot_type_counts: Path
    plot_missing_check: Path
    plot_missingness_heatmap: Path
    plot_log2fc_boxplot: Path
    plot_volcano_dawn: Path
    plot_volcano_dusk: Path
    plot_top_proteins: Path
    log_md: Path


def output_path(subdir: str, name: str, ext: str) -> Path:
    """Create a timestamped output path under outputs/."""
    ts = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    path = Path("outputs") / subdir / f"{ts}_{name}.{ext}"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def decode_annotation(value: Any) -> str:
    """Decode URL-encoded protein annotation text."""
    if pd.isna(value):
        return "unannotated"
    return unquote(str(value)).replace("_", " ")


def load_ptm_data() -> pd.DataFrame:
    """Load abundance-corrected PTM table and join annotation lookup."""
    ptm_path = Path("data/BRaVE_MED4_light-dark_dawn-dusk_all_ptms.tsv")
    ann_path = Path("data/protein_annotation_lookup.tsv")

    df = pd.read_csv(ptm_path, sep="\t")
    ann = pd.read_csv(ann_path, sep="\t")

    df["coord_key"] = df["UniProt"].astype(str).str.rsplit("_", n=1).str[0]
    df = df.merge(ann, left_on="coord_key", right_on="coordinates", how="left")
    df["protein_name"] = df["NCBI_annotation"].apply(decode_annotation)

    return df


def summarize_significance(df: pd.DataFrame, comparison: str) -> pd.DataFrame:
    """Summarize significant site counts by PTM type for a comparison."""
    adj_col = f"{comparison}_adj-p"
    if comparison not in df.columns or adj_col not in df.columns:
        return pd.DataFrame(columns=["Type", "direction", "count"])

    work = df[["Type", comparison, adj_col]].copy()
    work = work.replace([np.inf, -np.inf], np.nan).dropna()

    sig = work[work[adj_col] < 0.05].copy()
    # Fold-change columns can be ratio-scale or already log2-transformed.
    is_ratio_scale = bool(
        (work[comparison] > 0).all() and (work[comparison] >= 0).any()
    )
    if is_ratio_scale:
        sig["direction"] = np.where(sig[comparison] > 1.0, "increase", "decrease")
    else:
        sig["direction"] = np.where(sig[comparison] > 0.0, "increase", "decrease")

    out = (
        sig.groupby(["Type", "direction"], observed=True)
        .size()
        .reset_index(name="count")
        .sort_values(["Type", "direction"])
    )
    return out


def make_volcano(df: pd.DataFrame, comparison: str, title: str, out_path: Path) -> int:
    """Create a volcano-style plot and return significant point count."""
    adj_col = f"{comparison}_adj-p"
    plot_df = df[["Type", comparison, adj_col]].copy()
    plot_df = plot_df.replace([np.inf, -np.inf], np.nan).dropna()

    is_ratio_scale = bool((plot_df[comparison] > 0).all())
    if is_ratio_scale:
        plot_df["log2_fc"] = np.log2(plot_df[comparison])
    else:
        plot_df["log2_fc"] = plot_df[comparison]

    eps = np.finfo(float).tiny
    plot_df["neg_log10_adj_p"] = -np.log10(np.clip(plot_df[adj_col], eps, None))
    plot_df["significant"] = plot_df[adj_col] < 0.05

    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        data=plot_df,
        x="log2_fc",
        y="neg_log10_adj_p",
        hue="Type",
        style="significant",
        alpha=0.7,
        s=28,
        palette={"Ox": "#1f77b4", "Ac": "#ff7f0e"},
    )
    plt.axhline(-np.log10(0.05), color="black", linestyle="--", linewidth=1)
    plt.axvline(0, color="gray", linestyle=":", linewidth=1)
    plt.title(title)
    plt.xlabel("log2(Fold change)")
    plt.ylabel("-log10(adjusted p-value)")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()

    return int(plot_df["significant"].sum())


def main() -> None:
    """Run end-to-end EDA and save summary, plots, and log files."""
    request_text = (
        "Review copilot instructions and local skills, then perform full data exploration "
        "of abundance-corrected PTM data with concise summary and supporting plots."
    )
    plan_steps = [
        "Load abundance-corrected PTM table and join annotation lookup.",
        "Apply missingness-aware filter (missing_check == 0) for inferential summaries.",
        "Compute structure, quality, and significance summaries by PTM type and transition.",
        "Generate supporting visualization set.",
        "Write concise summary markdown and full execution log.",
    ]

    paths = ArtifactPaths(
        summary_md=output_path("artifacts", "abundance_corrected_eda_summary", "md"),
        filtered_tsv=output_path(
            "intermediate", "abundance_corrected_filtered_missingcheck0", "tsv"
        ),
        plot_type_counts=output_path(
            "artifacts", "plot_type_counts_prepost_filter", "png"
        ),
        plot_missing_check=output_path(
            "artifacts", "plot_missing_check_distribution", "png"
        ),
        plot_missingness_heatmap=output_path(
            "artifacts", "plot_condition_missingness_heatmap", "png"
        ),
        plot_log2fc_boxplot=output_path(
            "artifacts", "plot_log2fc_transition_boxplot", "png"
        ),
        plot_volcano_dawn=output_path(
            "artifacts", "plot_volcano_dawn_treatment", "png"
        ),
        plot_volcano_dusk=output_path(
            "artifacts", "plot_volcano_dusk_treatment", "png"
        ),
        plot_top_proteins=output_path(
            "artifacts", "plot_top_significant_proteins", "png"
        ),
        log_md=output_path("logs", "full_eda_abundance_corrected", "md"),
    )

    artifacts_for_log: list[str] = []
    errors: list[str] = []
    execution_notes: list[str] = []

    try:
        sns.set_theme(style="whitegrid", context="talk")

        df = load_ptm_data()
        execution_notes.append(
            f"Loaded abundance-corrected table with shape {df.shape[0]} rows x {df.shape[1]} columns."
        )

        # Required quality filter for confident statistics.
        df_f = df[df["missing_check"] == 0].copy()
        execution_notes.append(
            f"Applied missingness filter missing_check == 0, retaining {df_f.shape[0]} rows "
            f"({df_f.shape[0] / df.shape[0]:.1%} of total)."
        )

        df_f.to_csv(paths.filtered_tsv, sep="\t", index=False)
        artifacts_for_log.append(str(paths.filtered_tsv))

        base_counts = df["Type"].value_counts(dropna=False).rename("all_rows")
        filt_counts = df_f["Type"].value_counts(dropna=False).rename("filtered_rows")
        count_df = pd.concat([base_counts, filt_counts], axis=1).fillna(0).reset_index()
        count_df = count_df.rename(columns={"index": "Type"})
        count_df = count_df.melt(id_vars="Type", var_name="subset", value_name="count")

        plt.figure(figsize=(8, 5))
        sns.barplot(data=count_df, x="Type", y="count", hue="subset")
        plt.title("PTM type counts before and after missingness filter")
        plt.xlabel("PTM type")
        plt.ylabel("Site count")
        plt.tight_layout()
        plt.savefig(paths.plot_type_counts, dpi=300)
        plt.close()
        artifacts_for_log.append(str(paths.plot_type_counts))

        plt.figure(figsize=(8, 5))
        bins = np.arange(df["missing_check"].min(), df["missing_check"].max() + 2) - 0.5
        plt.hist(df["missing_check"], bins=bins, color="#4c78a8", edgecolor="white")
        plt.title("Distribution of missing_check across all PTM sites")
        plt.xlabel("missing_check (conditions with >50% replicate missingness)")
        plt.ylabel("Count")
        plt.tight_layout()
        plt.savefig(paths.plot_missing_check, dpi=300)
        plt.close()
        artifacts_for_log.append(str(paths.plot_missing_check))

        missingness_cols = [c for c in df.columns if c.endswith(" missingness")]
        miss_heat = (
            df_f.groupby("Type", observed=True)[missingness_cols]
            .mean()
            .T.rename_axis("condition")
        )
        plt.figure(figsize=(9, 4.8))
        sns.heatmap(
            miss_heat,
            cmap="mako",
            annot=True,
            fmt=".2f",
            cbar_kws={"label": "mean missing replicates"},
        )
        plt.title("Mean per-condition missingness after filtering")
        plt.xlabel("PTM type")
        plt.ylabel("Condition")
        plt.tight_layout()
        plt.savefig(paths.plot_missingness_heatmap, dpi=300)
        plt.close()
        artifacts_for_log.append(str(paths.plot_missingness_heatmap))

        transitions = {
            "Dawn treatment (light pulse)": "Light_0630/Dark_0600",
            "Dawn control (dark continuation)": "Dark_0630/Dark_0600",
            "Dusk treatment (dark pulse)": "Dark_1830/Light_1800",
            "Dusk control (light continuation)": "Light_1830/Light_1800",
        }

        long_fc_rows = []
        for label, col in transitions.items():
            if col in df_f.columns:
                sub = df_f[["Type", col]].copy()
                sub["transition"] = label
                sub = sub.replace([np.inf, -np.inf], np.nan).dropna()
                is_ratio_scale = bool((sub[col] > 0).all())
                if is_ratio_scale:
                    sub["log2_fc"] = np.log2(sub[col])
                else:
                    sub["log2_fc"] = sub[col]
                long_fc_rows.append(sub[["Type", "transition", "log2_fc"]])

        long_fc = pd.concat(long_fc_rows, ignore_index=True)
        plt.figure(figsize=(11, 5.5))
        sns.boxplot(
            data=long_fc, x="transition", y="log2_fc", hue="Type", showfliers=False
        )
        plt.axhline(0, color="gray", linestyle=":", linewidth=1)
        plt.title("log2 fold-change distributions by transition and PTM type")
        plt.xlabel("Transition")
        plt.ylabel("log2(fold change)")
        plt.xticks(rotation=15, ha="right")
        plt.tight_layout()
        plt.savefig(paths.plot_log2fc_boxplot, dpi=300)
        plt.close()
        artifacts_for_log.append(str(paths.plot_log2fc_boxplot))

        dawn_sig_n = make_volcano(
            df_f,
            "Light_0630/Dark_0600",
            "Dawn treatment volcano: Light_0630 vs Dark_0600",
            paths.plot_volcano_dawn,
        )
        dusk_sig_n = make_volcano(
            df_f,
            "Dark_1830/Light_1800",
            "Dusk treatment volcano: Dark_1830 vs Light_1800",
            paths.plot_volcano_dusk,
        )
        artifacts_for_log.append(str(paths.plot_volcano_dawn))
        artifacts_for_log.append(str(paths.plot_volcano_dusk))

        # Top proteins by significant sites across two treatment comparisons.
        sig_mask = (df_f["Light_0630/Dark_0600_adj-p"] < 0.05) | (
            df_f["Dark_1830/Light_1800_adj-p"] < 0.05
        )
        top_df = (
            df_f.loc[sig_mask, ["protein_name"]]
            .value_counts()
            .reset_index(name="significant_sites")
            .head(12)
        )

        plt.figure(figsize=(9, 6))
        sns.barplot(
            data=top_df, x="significant_sites", y="protein_name", color="#2a9d8f"
        )
        plt.title("Top annotated proteins by significant PTM site count")
        plt.xlabel("Significant sites (adj-p < 0.05 in dawn or dusk treatment)")
        plt.ylabel("Protein annotation")
        plt.tight_layout()
        plt.savefig(paths.plot_top_proteins, dpi=300)
        plt.close()
        artifacts_for_log.append(str(paths.plot_top_proteins))

        # Numerical summaries.
        unique_proteins_all = df["UniProt"].nunique(dropna=True)
        unique_proteins_filtered = df_f["UniProt"].nunique(dropna=True)
        annotation_rate = float((df["protein_name"] != "unannotated").mean())

        sig_tables = {
            "Dawn treatment": summarize_significance(df_f, "Light_0630/Dark_0600"),
            "Dusk treatment": summarize_significance(df_f, "Dark_1830/Light_1800"),
            "Dawn endpoint (Light vs Dark at 06:30)": summarize_significance(
                df_f, "Light_0630/Dark_0630"
            ),
            "Dusk endpoint (Dark vs Light at 18:30)": summarize_significance(
                df_f, "Dark_1830/Light_1830"
            ),
        }

        summary_lines = [
            "# Abundance-Corrected PTM Data Exploration (Concise Summary)",
            "",
            f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            "",
            "## Dataset Snapshot",
            f"- Input file: data/BRaVE_MED4_light-dark_dawn-dusk_all_ptms.tsv",
            f"- Rows x columns: {df.shape[0]:,} x {df.shape[1]:,}",
            f"- Unique proteins (all rows): {unique_proteins_all:,}",
            f"- Unique proteins (filtered rows): {unique_proteins_filtered:,}",
            f"- Annotation coverage (non-unannotated rows): {annotation_rate:.1%}",
            "",
            "## Filtering and Data Quality",
            "- Statistical summaries were based on missing_check == 0 (required quality filter).",
            f"- Retained rows after filter: {df_f.shape[0]:,} ({df_f.shape[0] / df.shape[0]:.1%} of all rows)",
            f"- Mean total missingness (all rows): {df['Total missingness'].mean():.2f}",
            f"- Mean total missingness (filtered rows): {df_f['Total missingness'].mean():.2f}",
            "",
            "## Key Transition-Level Signals (adj-p < 0.05)",
            f"- Dawn treatment (Light_0630/Dark_0600): {dawn_sig_n:,} significant sites",
            f"- Dusk treatment (Dark_1830/Light_1800): {dusk_sig_n:,} significant sites",
            "",
            "## Significant Site Breakdown by PTM Type",
        ]

        for title, tab in sig_tables.items():
            summary_lines.append(f"### {title}")
            if tab.empty:
                summary_lines.append("- No results (required columns not present).")
            else:
                for _, row in tab.iterrows():
                    summary_lines.append(
                        f"- {row['Type']} {row['direction']}: {int(row['count']):,} sites"
                    )
            summary_lines.append("")

        summary_lines.extend(
            [
                "## Supporting Plot Artifacts",
                f"- {paths.plot_type_counts}",
                f"- {paths.plot_missing_check}",
                f"- {paths.plot_missingness_heatmap}",
                f"- {paths.plot_log2fc_boxplot}",
                f"- {paths.plot_volcano_dawn}",
                f"- {paths.plot_volcano_dusk}",
                f"- {paths.plot_top_proteins}",
            ]
        )

        paths.summary_md.write_text("\n".join(summary_lines), encoding="utf-8")
        artifacts_for_log.append(str(paths.summary_md))

        execution_notes.append(
            "Generated summary markdown and 7 supporting plots spanning structure, quality, and transition-specific PTM responses."
        )

    except Exception as exc:  # pragma: no cover
        errors.append(f"{type(exc).__name__}: {exc}")
        errors.append(traceback.format_exc())

    log_lines = [
        "Task: Full EDA of abundance-corrected MED4 PTM data",
        "",
        f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        f"Request: {request_text}",
        "",
        "Plan",
        "",
    ]

    for i, step in enumerate(plan_steps, start=1):
        log_lines.append(f"Step {i}: {step}")

    log_lines.extend(["", "Execution", ""])

    for i, note in enumerate(execution_notes, start=1):
        log_lines.extend(
            [
                f"Step {i}: Analysis stage",
                "",
                "What was done",
                note,
                "",
                "Result or observation",
                note,
                "",
            ]
        )

    log_lines.extend(["Results", ""])

    if errors:
        log_lines.append("Run completed with errors; see Errors / Warnings.")
    else:
        log_lines.append("EDA completed successfully with summary and plots generated.")

    log_lines.extend(["", "Artifacts Produced", ""])

    if artifacts_for_log:
        log_lines.extend(artifacts_for_log)
    else:
        log_lines.append("None")

    log_lines.extend(["", "Errors / Warnings", ""])

    if errors:
        log_lines.extend(errors)
    else:
        log_lines.append("None")

    paths.log_md.write_text("\n".join(log_lines), encoding="utf-8")

    if errors:
        raise RuntimeError(f"EDA failed. See log: {paths.log_md}")

    print("EDA complete.")
    print(f"Summary: {paths.summary_md}")
    print(f"Log: {paths.log_md}")
    print("Plots:")
    print(paths.plot_type_counts)
    print(paths.plot_missing_check)
    print(paths.plot_missingness_heatmap)
    print(paths.plot_log2fc_boxplot)
    print(paths.plot_volcano_dawn)
    print(paths.plot_volcano_dusk)
    print(paths.plot_top_proteins)


if __name__ == "__main__":
    main()
