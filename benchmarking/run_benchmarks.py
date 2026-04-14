#!/usr/bin/env python
"""
Main benchmark runner for spaceprime.

Runs timing and memory benchmarks across various parameter combinations.

Usage:
    python run_benchmarks.py --config quick   # Fast test run
    python run_benchmarks.py --config default # Full suite
    python run_benchmarks.py --output results/benchmark_2024.csv
"""

import argparse
import csv
import os
import sys
import time
import uuid
from datetime import datetime
from itertools import product
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

import msprime

from spaceprime import spDemography
from spaceprime.utilities import calc_migration_matrix

from benchmark_utils import (
    MemoryTracker,
    create_sample_dict,
    generate_synthetic_demes,
    time_function,
)
from config import get_config


# CSV column headers
CSV_COLUMNS = [
    "run_id",
    "timestamp",
    "grid_rows",
    "grid_cols",
    "num_demes",
    "max_local_size",
    "migration_rate",
    "time_slices",
    "timesteps",
    "replicate",
    "migration_matrix_time_sec",
    "demography_setup_time_sec",
    "simulation_time_sec",
    "total_time_sec",
    "peak_memory_mb",
    "num_samples",
    "completed",
    "error",
]


def run_single_benchmark(
    grid_shape: tuple,
    max_local_size: int,
    migration_rate: float,
    time_slices: int,
    timesteps: int,
    replicate: int,
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Run a single benchmark for a given parameter combination.

    Parameters
    ----------
    grid_shape : tuple
        Shape of the grid (rows, cols).
    max_local_size : int
        Maximum deme size.
    migration_rate : float
        Migration rate between demes.
    time_slices : int
        Number of time slices.
    timesteps : int
        Timesteps between demographic events.
    replicate : int
        Replicate number.
    verbose : bool, optional
        Print progress. Default is True.

    Returns
    -------
    Dict[str, Any]
        Dictionary containing benchmark results.
    """
    run_id = str(uuid.uuid4())[:8]
    timestamp = datetime.now().isoformat()
    rows, cols = grid_shape
    num_demes = rows * cols

    result = {
        "run_id": run_id,
        "timestamp": timestamp,
        "grid_rows": rows,
        "grid_cols": cols,
        "num_demes": num_demes,
        "max_local_size": max_local_size,
        "migration_rate": migration_rate,
        "time_slices": time_slices,
        "timesteps": timesteps,
        "replicate": replicate,
        "migration_matrix_time_sec": None,
        "demography_setup_time_sec": None,
        "simulation_time_sec": None,
        "total_time_sec": None,
        "peak_memory_mb": None,
        "num_samples": None,
        "completed": False,
        "error": None,
    }

    if verbose:
        print(
            f"  Running: {rows}x{cols} grid, size={max_local_size}, "
            f"mig={migration_rate}, slices={time_slices}, rep={replicate}"
        )

    try:
        total_start = time.time()

        with MemoryTracker() as tracker:
            # Generate synthetic demes
            demes = generate_synthetic_demes(
                grid_shape=grid_shape,
                max_size=max_local_size,
                time_slices=time_slices,
                pattern="random",
                seed=42 + replicate,
            )

            # Get 2D version for migration matrix
            if len(demes.shape) == 3:
                demes_2d = demes[0]
            else:
                demes_2d = demes

            # Time calc_migration_matrix
            _, mig_time = time_function(
                calc_migration_matrix, demes_2d, migration_rate, scale=True
            )
            result["migration_matrix_time_sec"] = mig_time

            # Time spDemography.stepping_stone_2d
            demo = spDemography()
            demo_start = time.time()
            demo.stepping_stone_2d(
                d=demes,
                rate=migration_rate,
                scale=True,
                timesteps=timesteps if time_slices > 1 else None,
            )
            demo.add_ancestral_populations(anc_sizes=[10000], merge_time=10000)
            demo_time = time.time() - demo_start
            result["demography_setup_time_sec"] = demo_time

            # Create sample dict and run simulation
            sample_dict = create_sample_dict(demes, samples_per_deme=2, min_deme_size=10)
            result["num_samples"] = sum(sample_dict.values())

            # Only run simulation if we have samples
            if sample_dict:
                sim_start = time.time()
                ts = msprime.sim_ancestry(
                    samples=sample_dict,
                    demography=demo,
                    sequence_length=10000,  # Short sequence for benchmarking
                    recombination_rate=0,
                    ploidy=2,
                    random_seed=42 + replicate,
                    record_provenance=False,
                )
                sim_time = time.time() - sim_start
                result["simulation_time_sec"] = sim_time
            else:
                result["simulation_time_sec"] = 0.0
                result["error"] = "No valid demes to sample from"

        result["peak_memory_mb"] = tracker.peak_mb
        result["total_time_sec"] = time.time() - total_start
        result["completed"] = True

    except Exception as e:
        result["error"] = str(e)
        result["completed"] = False
        if verbose:
            print(f"    Error: {e}")

    return result


def run_benchmarks(
    config_name: str = "default",
    output_path: Optional[str] = None,
    verbose: bool = True,
) -> List[Dict[str, Any]]:
    """
    Run the full benchmark suite.

    Parameters
    ----------
    config_name : str, optional
        Configuration name. Options: 'quick', 'default', 'full'. Default is 'default'.
    output_path : str, optional
        Path to save results CSV. Default is None (auto-generated).
    verbose : bool, optional
        Print progress. Default is True.

    Returns
    -------
    List[Dict[str, Any]]
        List of benchmark results.
    """
    config = get_config(config_name)

    if output_path is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_dir = Path(__file__).parent / "results"
        output_dir.mkdir(exist_ok=True)
        output_path = output_dir / f"benchmark_{config_name}_{timestamp}.csv"
    else:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

    # Calculate total number of runs
    total_runs = (
        len(config["deme_sizes"])
        * len(config["grid_dimensions"])
        * len(config["migration_rates"])
        * len(config["time_slices"])
        * len(config["timesteps"])
        * config["replicates"]
    )

    if verbose:
        print(f"Starting benchmark suite: {config_name}")
        print(f"Total parameter combinations: {total_runs}")
        print(f"Output: {output_path}")
        print("-" * 60)

    results = []
    run_count = 0

    # Create/open CSV file
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()

        # Iterate through all parameter combinations
        for grid_dim in config["grid_dimensions"]:
            for max_size in config["deme_sizes"]:
                for mig_rate in config["migration_rates"]:
                    for time_slices in config["time_slices"]:
                        for timesteps_val in config["timesteps"]:
                            for rep in range(config["replicates"]):
                                run_count += 1

                                if verbose:
                                    print(f"[{run_count}/{total_runs}]", end="")

                                result = run_single_benchmark(
                                    grid_shape=grid_dim,
                                    max_local_size=max_size,
                                    migration_rate=mig_rate,
                                    time_slices=time_slices,
                                    timesteps=timesteps_val,
                                    replicate=rep,
                                    verbose=verbose,
                                )

                                results.append(result)
                                writer.writerow(result)
                                f.flush()  # Ensure data is written incrementally

                                if verbose and result["completed"]:
                                    print(
                                        f"    Total: {result['total_time_sec']:.2f}s, "
                                        f"Memory: {result['peak_memory_mb']:.1f}MB"
                                    )

    if verbose:
        print("-" * 60)
        print(f"Benchmark complete. Results saved to: {output_path}")
        completed = sum(1 for r in results if r["completed"])
        print(f"Completed: {completed}/{len(results)}")

    return results


def main():
    """Main entry point for benchmark runner."""
    parser = argparse.ArgumentParser(
        description="Run spaceprime benchmarks",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_benchmarks.py --config quick
  python run_benchmarks.py --config default --output results/my_benchmark.csv
  python run_benchmarks.py --config full --quiet
        """,
    )

    parser.add_argument(
        "--config",
        "-c",
        type=str,
        default="default",
        choices=["quick", "default", "full"],
        help="Benchmark configuration (default: default)",
    )

    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default=None,
        help="Output CSV path (default: auto-generated in results/)",
    )

    parser.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="Suppress progress output",
    )

    args = parser.parse_args()

    run_benchmarks(
        config_name=args.config,
        output_path=args.output,
        verbose=not args.quiet,
    )


if __name__ == "__main__":
    main()
