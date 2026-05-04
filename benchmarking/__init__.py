"""
spaceprime benchmarking module.

Tools for measuring performance of spaceprime simulations.
"""

from .config import get_config, DEFAULT_CONFIG, QUICK_CONFIG, FULL_CONFIG
from .benchmark_utils import (
    generate_synthetic_demes,
    create_sample_dict,
    timed,
    timer,
    MemoryTracker,
    time_function,
)

__all__ = [
    "get_config",
    "DEFAULT_CONFIG",
    "QUICK_CONFIG",
    "FULL_CONFIG",
    "generate_synthetic_demes",
    "create_sample_dict",
    "timed",
    "timer",
    "MemoryTracker",
    "time_function",
]
