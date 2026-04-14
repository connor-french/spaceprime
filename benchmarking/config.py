"""
Configuration settings for spaceprime benchmarking.

Defines parameter ranges for different benchmark configurations.
"""

from typing import Dict, List, Tuple, Any

# Default configuration for comprehensive benchmarking
DEFAULT_CONFIG: Dict[str, Any] = {
    "deme_sizes": [100, 500, 1000, 5000, 10000],
    "grid_dimensions": [(5, 5), (10, 10), (20, 20), (30, 30)],
    "migration_rates": [0.001, 0.01, 0.1],
    "time_slices": [1, 5, 10, 20],
    "timesteps": [100, 1000],
    "replicates": 3,
}

# Quick configuration for testing the benchmark setup
QUICK_CONFIG: Dict[str, Any] = {
    "deme_sizes": [1000],
    "grid_dimensions": [(5, 5), (10, 10)],
    "migration_rates": [0.01],
    "time_slices": [1, 5],
    "timesteps": [100],
    "replicates": 1,
}

# Full configuration for exhaustive benchmarking
FULL_CONFIG: Dict[str, Any] = {
    "deme_sizes": [100, 500, 1000, 2500, 5000, 7500, 10000],
    "grid_dimensions": [(5, 5), (10, 10), (15, 15), (20, 20), (25, 25), (30, 30)],
    "migration_rates": [0.0001, 0.001, 0.01, 0.05, 0.1],
    "time_slices": [1, 2, 5, 10, 15, 20],
    "timesteps": [100, 500, 1000],
    "replicates": 5,
}


def get_config(config_name: str = "default") -> Dict[str, Any]:
    """
    Get benchmark configuration by name.

    Parameters
    ----------
    config_name : str
        Name of the configuration. Options: 'quick', 'default', 'full'.

    Returns
    -------
    Dict[str, Any]
        Configuration dictionary with parameter ranges.

    Raises
    ------
    ValueError
        If config_name is not recognized.
    """
    configs = {
        "quick": QUICK_CONFIG,
        "default": DEFAULT_CONFIG,
        "full": FULL_CONFIG,
    }

    if config_name not in configs:
        raise ValueError(
            f"Unknown config '{config_name}'. Options: {list(configs.keys())}"
        )

    return configs[config_name]
