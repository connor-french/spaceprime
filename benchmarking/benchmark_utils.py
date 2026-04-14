"""
Utility functions for spaceprime benchmarking.

Provides synthetic data generation, timing utilities, and memory tracking.
"""

import time
import tracemalloc
from contextlib import contextmanager
from functools import wraps
from typing import Callable, Dict, List, Optional, Tuple, Any

import numpy as np


def generate_synthetic_demes(
    grid_shape: Tuple[int, int],
    max_size: int,
    time_slices: int = 1,
    pattern: str = "uniform",
    seed: Optional[int] = None,
) -> np.ndarray:
    """
    Generate synthetic deme arrays for benchmarking.

    Parameters
    ----------
    grid_shape : Tuple[int, int]
        Shape of the grid (rows, cols).
    max_size : int
        Maximum deme size.
    time_slices : int, optional
        Number of time slices. If > 1, returns a 3D array. Default is 1.
    pattern : str, optional
        Pattern for generating deme sizes. Options:
        - "uniform": All demes have max_size
        - "random": Random sizes between 0 and max_size
        - "gradient": Linear gradient across the grid
        Default is "uniform".
    seed : int, optional
        Random seed for reproducibility. Default is None.

    Returns
    -------
    np.ndarray
        2D array (rows, cols) if time_slices=1, else 3D array (time, rows, cols).
    """
    if seed is not None:
        np.random.seed(seed)

    rows, cols = grid_shape

    if pattern == "uniform":
        base = np.full((rows, cols), max_size, dtype=float)
    elif pattern == "random":
        # Random values between 10% and 100% of max_size
        base = np.random.uniform(0.1, 1.0, (rows, cols)) * max_size
        base = np.ceil(base)
    elif pattern == "gradient":
        # Linear gradient from 10% to 100% across the grid
        x = np.linspace(0.1, 1.0, cols)
        y = np.linspace(0.1, 1.0, rows)
        xx, yy = np.meshgrid(x, y)
        base = ((xx + yy) / 2) * max_size
        base = np.ceil(base)
    else:
        raise ValueError(f"Unknown pattern: {pattern}. Use 'uniform', 'random', or 'gradient'")

    # Set some cells to near-zero to simulate unsuitable habitat
    mask = np.random.random((rows, cols)) > 0.85
    base[mask] = 1e-10

    if time_slices == 1:
        return base

    # For multiple time slices, create variations
    result = np.zeros((time_slices, rows, cols))
    result[0] = base

    for t in range(1, time_slices):
        # Each time slice is a modified version (simulating temporal changes)
        variation = np.random.uniform(0.8, 1.2, (rows, cols))
        result[t] = np.clip(base * variation, 1e-10, max_size)
        result[t] = np.ceil(result[t])
        # Preserve the unsuitable habitat mask
        result[t][mask] = 1e-10

    return result


def create_sample_dict(
    demes: np.ndarray,
    samples_per_deme: int = 2,
    min_deme_size: int = 10,
) -> Dict[str, int]:
    """
    Create a sample dictionary matching coords_to_sample_dict output format.

    Parameters
    ----------
    demes : np.ndarray
        2D or 3D array of deme sizes. If 3D, uses the first time slice.
    samples_per_deme : int, optional
        Number of samples per deme. Default is 2.
    min_deme_size : int, optional
        Minimum deme size to sample from. Default is 10.

    Returns
    -------
    Dict[str, int]
        Dictionary mapping deme names to sample counts.
    """
    if len(demes.shape) == 3:
        demes_2d = demes[0]
    else:
        demes_2d = demes

    sample_dict = {}
    rows, cols = demes_2d.shape

    for i in range(rows):
        for j in range(cols):
            if demes_2d[i, j] >= min_deme_size:
                deme_name = f"deme_{i}_{j}"
                sample_dict[deme_name] = samples_per_deme

    return sample_dict


def timed(func: Callable) -> Callable:
    """
    Decorator to measure execution time of a function.

    Returns a tuple of (result, elapsed_time_seconds).
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        elapsed_time = time.time() - start_time
        return result, elapsed_time
    return wrapper


@contextmanager
def timer():
    """
    Context manager for timing code blocks.

    Yields a dictionary that will contain 'elapsed' key after the block completes.

    Example
    -------
    >>> with timer() as t:
    ...     do_something()
    >>> print(f"Elapsed: {t['elapsed']:.2f}s")
    """
    timing = {}
    start = time.time()
    try:
        yield timing
    finally:
        timing['elapsed'] = time.time() - start


class MemoryTracker:
    """
    Context manager for tracking memory usage using tracemalloc.

    Attributes
    ----------
    peak_mb : float
        Peak memory usage in megabytes after exiting the context.
    current_mb : float
        Current memory usage in megabytes after exiting the context.

    Example
    -------
    >>> with MemoryTracker() as tracker:
    ...     # allocate memory
    ...     data = [i for i in range(1000000)]
    >>> print(f"Peak: {tracker.peak_mb:.2f} MB")
    """

    def __init__(self):
        self.peak_mb = 0.0
        self.current_mb = 0.0

    def __enter__(self):
        tracemalloc.start()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        self.peak_mb = peak / (1024 * 1024)
        self.current_mb = current / (1024 * 1024)
        return False


def time_function(func: Callable, *args, **kwargs) -> Tuple[Any, float]:
    """
    Time a function call.

    Parameters
    ----------
    func : Callable
        Function to time.
    *args : Any
        Positional arguments to pass to the function.
    **kwargs : Any
        Keyword arguments to pass to the function.

    Returns
    -------
    Tuple[Any, float]
        Tuple of (result, elapsed_time_seconds).
    """
    start = time.time()
    result = func(*args, **kwargs)
    elapsed = time.time() - start
    return result, elapsed


def estimate_memory_usage(grid_shape: Tuple[int, int], time_slices: int = 1) -> float:
    """
    Estimate memory usage for a given grid configuration.

    Parameters
    ----------
    grid_shape : Tuple[int, int]
        Shape of the grid (rows, cols).
    time_slices : int, optional
        Number of time slices. Default is 1.

    Returns
    -------
    float
        Estimated memory usage in MB.
    """
    rows, cols = grid_shape
    n_demes = rows * cols

    # Migration matrix: n_demes x n_demes, float64
    migration_matrix_mb = (n_demes ** 2) * 8 / (1024 * 1024)

    # Demes array
    demes_mb = (rows * cols * time_slices) * 8 / (1024 * 1024)

    # Rough estimate of msprime overhead (populations, events)
    msprime_overhead_mb = n_demes * 0.001  # ~1KB per population

    return migration_matrix_mb + demes_mb + msprime_overhead_mb
