#!/usr/bin/env python
"""
Analyze and visualize spaceprime benchmark results.

This script generates summary statistics and visualizations from benchmark CSV files.

Usage:
    python analyze_results.py results/benchmark_default_20240413.csv
    python analyze_results.py results/*.csv --output figures/
"""

import argparse
from pathlib import Path
from typing import List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def load_results(file_paths: List[str]) -> pd.DataFrame:
    """
    Load and combine benchmark results from one or more CSV files.

    Parameters
    ----------
    file_paths : List[str]
        List of paths to benchmark CSV files.

    Returns
    -------
    pd.DataFrame
        Combined benchmark results.
    """
    dfs = []
    for path in file_paths:
        df = pd.read_csv(path)
        df["source_file"] = Path(path).name
        dfs.append(df)

    combined = pd.concat(dfs, ignore_index=True)

    # Filter to completed runs
    combined = combined[combined["completed"] == True].copy()

    return combined


def summarize_results(df: pd.DataFrame) -> pd.DataFrame:
    """
    Generate summary statistics for benchmark results.

    Parameters
    ----------
    df : pd.DataFrame
        Benchmark results DataFrame.

    Returns
    -------
    pd.DataFrame
        Summary statistics grouped by parameter combinations.
    """
    group_cols = [
        "grid_rows",
        "grid_cols",
        "num_demes",
        "max_local_size",
        "migration_rate",
        "time_slices",
        "timesteps",
    ]

    summary = (
        df.groupby(group_cols)
        .agg(
            {
                "migration_matrix_time_sec": ["mean", "std"],
                "demography_setup_time_sec": ["mean", "std"],
                "simulation_time_sec": ["mean", "std"],
                "total_time_sec": ["mean", "std"],
                "peak_memory_mb": ["mean", "std", "max"],
                "num_samples": "first",
                "replicate": "count",
            }
        )
        .reset_index()
    )

    # Flatten column names
    summary.columns = [
        "_".join(col).strip("_") if isinstance(col, tuple) else col
        for col in summary.columns
    ]

    return summary


def plot_scaling_by_demes(df: pd.DataFrame, output_dir: Optional[Path] = None) -> plt.Figure:
    """
    Plot log-log scaling of setup time vs number of demes.

    Parameters
    ----------
    df : pd.DataFrame
        Benchmark results DataFrame.
    output_dir : Path, optional
        Directory to save figure. Default is None (don't save).

    Returns
    -------
    plt.Figure
        Matplotlib figure.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Group by num_demes
    grouped = df.groupby("num_demes").agg(
        {
            "demography_setup_time_sec": ["mean", "std"],
            "total_time_sec": ["mean", "std"],
        }
    )

    x = grouped.index.values

    # Left plot: Demography setup time
    ax = axes[0]
    y = grouped["demography_setup_time_sec"]["mean"].values
    yerr = grouped["demography_setup_time_sec"]["std"].values
    ax.errorbar(x, y, yerr=yerr, fmt="o-", capsize=3, label="Setup time")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of Demes")
    ax.set_ylabel("Time (seconds)")
    ax.set_title("Demography Setup Time vs. Grid Size")
    ax.grid(True, alpha=0.3)

    # Right plot: Total time
    ax = axes[1]
    y = grouped["total_time_sec"]["mean"].values
    yerr = grouped["total_time_sec"]["std"].values
    ax.errorbar(x, y, yerr=yerr, fmt="s-", capsize=3, color="C1", label="Total time")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of Demes")
    ax.set_ylabel("Time (seconds)")
    ax.set_title("Total Benchmark Time vs. Grid Size")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_dir:
        fig.savefig(output_dir / "scaling_by_demes.png", dpi=150, bbox_inches="tight")

    return fig


def plot_heatmap_grid_time(df: pd.DataFrame, output_dir: Optional[Path] = None) -> plt.Figure:
    """
    Plot heatmap of grid size × time slices → total time.

    Parameters
    ----------
    df : pd.DataFrame
        Benchmark results DataFrame.
    output_dir : Path, optional
        Directory to save figure. Default is None (don't save).

    Returns
    -------
    plt.Figure
        Matplotlib figure.
    """
    # Create pivot table: grid size string vs time slices
    df = df.copy()
    df["grid_size"] = df["grid_rows"].astype(str) + "x" + df["grid_cols"].astype(str)

    pivot = df.pivot_table(
        values="total_time_sec",
        index="grid_size",
        columns="time_slices",
        aggfunc="mean",
    )

    fig, ax = plt.subplots(figsize=(10, 6))

    im = ax.imshow(pivot.values, aspect="auto", cmap="YlOrRd")
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns)
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    ax.set_xlabel("Time Slices")
    ax.set_ylabel("Grid Size")
    ax.set_title("Mean Total Time (seconds)")

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Time (s)")

    # Add value annotations
    for i in range(len(pivot.index)):
        for j in range(len(pivot.columns)):
            value = pivot.values[i, j]
            if not np.isnan(value):
                text_color = "white" if value > pivot.values.max() / 2 else "black"
                ax.text(j, i, f"{value:.2f}", ha="center", va="center", color=text_color)

    plt.tight_layout()

    if output_dir:
        fig.savefig(output_dir / "heatmap_grid_time.png", dpi=150, bbox_inches="tight")

    return fig


def plot_parameter_effects(df: pd.DataFrame, output_dir: Optional[Path] = None) -> plt.Figure:
    """
    Plot effect of individual parameters on timing.

    Parameters
    ----------
    df : pd.DataFrame
        Benchmark results DataFrame.
    output_dir : Path, optional
        Directory to save figure. Default is None (don't save).

    Returns
    -------
    plt.Figure
        Matplotlib figure.
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Plot 1: Effect of max_local_size
    ax = axes[0, 0]
    grouped = df.groupby("max_local_size")["total_time_sec"].agg(["mean", "std"])
    ax.errorbar(
        grouped.index,
        grouped["mean"],
        yerr=grouped["std"],
        fmt="o-",
        capsize=3,
    )
    ax.set_xlabel("Max Local Size")
    ax.set_ylabel("Total Time (s)")
    ax.set_title("Effect of Deme Size")
    ax.grid(True, alpha=0.3)

    # Plot 2: Effect of migration rate
    ax = axes[0, 1]
    grouped = df.groupby("migration_rate")["total_time_sec"].agg(["mean", "std"])
    ax.errorbar(
        grouped.index,
        grouped["mean"],
        yerr=grouped["std"],
        fmt="s-",
        capsize=3,
        color="C1",
    )
    ax.set_xscale("log")
    ax.set_xlabel("Migration Rate")
    ax.set_ylabel("Total Time (s)")
    ax.set_title("Effect of Migration Rate")
    ax.grid(True, alpha=0.3)

    # Plot 3: Effect of time slices
    ax = axes[1, 0]
    grouped = df.groupby("time_slices")["total_time_sec"].agg(["mean", "std"])
    ax.errorbar(
        grouped.index,
        grouped["mean"],
        yerr=grouped["std"],
        fmt="^-",
        capsize=3,
        color="C2",
    )
    ax.set_xlabel("Time Slices")
    ax.set_ylabel("Total Time (s)")
    ax.set_title("Effect of Time Slices")
    ax.grid(True, alpha=0.3)

    # Plot 4: Memory scaling
    ax = axes[1, 1]
    grouped = df.groupby("num_demes")["peak_memory_mb"].agg(["mean", "std", "max"])
    ax.errorbar(
        grouped.index,
        grouped["mean"],
        yerr=grouped["std"],
        fmt="d-",
        capsize=3,
        color="C3",
        label="Mean",
    )
    ax.plot(grouped.index, grouped["max"], "x--", color="C4", label="Max")
    ax.set_xscale("log")
    ax.set_xlabel("Number of Demes")
    ax.set_ylabel("Peak Memory (MB)")
    ax.set_title("Memory Scaling")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_dir:
        fig.savefig(output_dir / "parameter_effects.png", dpi=150, bbox_inches="tight")

    return fig


def plot_time_breakdown(df: pd.DataFrame, output_dir: Optional[Path] = None) -> plt.Figure:
    """
    Plot stacked bar chart showing time breakdown by component.

    Parameters
    ----------
    df : pd.DataFrame
        Benchmark results DataFrame.
    output_dir : Path, optional
        Directory to save figure. Default is None (don't save).

    Returns
    -------
    plt.Figure
        Matplotlib figure.
    """
    fig, ax = plt.subplots(figsize=(12, 6))

    df = df.copy()
    df["grid_size"] = df["grid_rows"].astype(str) + "x" + df["grid_cols"].astype(str)

    # Group by grid size
    grouped = df.groupby("grid_size").agg(
        {
            "migration_matrix_time_sec": "mean",
            "demography_setup_time_sec": "mean",
            "simulation_time_sec": "mean",
        }
    )

    x = range(len(grouped))
    width = 0.8

    ax.bar(
        x,
        grouped["migration_matrix_time_sec"],
        width,
        label="Migration Matrix",
    )
    ax.bar(
        x,
        grouped["demography_setup_time_sec"],
        width,
        bottom=grouped["migration_matrix_time_sec"],
        label="Demography Setup",
    )
    ax.bar(
        x,
        grouped["simulation_time_sec"],
        width,
        bottom=grouped["migration_matrix_time_sec"] + grouped["demography_setup_time_sec"],
        label="Simulation",
    )

    ax.set_xticks(x)
    ax.set_xticklabels(grouped.index, rotation=45, ha="right")
    ax.set_xlabel("Grid Size")
    ax.set_ylabel("Time (seconds)")
    ax.set_title("Time Breakdown by Component")
    ax.legend()

    plt.tight_layout()

    if output_dir:
        fig.savefig(output_dir / "time_breakdown.png", dpi=150, bbox_inches="tight")

    return fig


def generate_report(df: pd.DataFrame, output_dir: Path) -> None:
    """
    Generate a full analysis report with all plots and summary statistics.

    Parameters
    ----------
    df : pd.DataFrame
        Benchmark results DataFrame.
    output_dir : Path
        Directory to save report files.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Generate summary statistics
    summary = summarize_results(df)
    summary.to_csv(output_dir / "summary_statistics.csv", index=False)

    # Generate plots
    plot_scaling_by_demes(df, output_dir)
    plot_heatmap_grid_time(df, output_dir)
    plot_parameter_effects(df, output_dir)
    plot_time_breakdown(df, output_dir)

    print(f"Report generated in: {output_dir}")
    print(f"  - summary_statistics.csv")
    print(f"  - scaling_by_demes.png")
    print(f"  - heatmap_grid_time.png")
    print(f"  - parameter_effects.png")
    print(f"  - time_breakdown.png")


def main():
    """Main entry point for analysis script."""
    parser = argparse.ArgumentParser(
        description="Analyze spaceprime benchmark results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python analyze_results.py results/benchmark_default_20240413.csv
  python analyze_results.py results/*.csv --output figures/
        """,
    )

    parser.add_argument(
        "files",
        nargs="+",
        help="Path(s) to benchmark CSV file(s)",
    )

    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default="figures",
        help="Output directory for figures (default: figures/)",
    )

    parser.add_argument(
        "--show",
        action="store_true",
        help="Display plots interactively",
    )

    args = parser.parse_args()

    # Load results
    df = load_results(args.files)
    print(f"Loaded {len(df)} benchmark results from {len(args.files)} file(s)")

    # Generate report
    output_dir = Path(args.output)
    generate_report(df, output_dir)

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
