"""
spaceprime is a python package to facilitate spatially explicit coalescent modeling in msprime.
Please see the documentation at https://connor-french.github.io/spaceprime
for more information.
"""

__author__ = """Connor French"""
__email__ = "french.connor.m@gmail.com"
__version__ = "0.0.6"

__all__ = [
    "__version__",
    # utilities
    "create_raster",
    "raster_to_demes",
    "calc_migration_matrix",
    "split_landscape_by_pop",
    "mtp_thresh_from_coords",
    "coords_to_sample_dict",
    "anc_to_deme_dict",
    "coords_to_deme_dict",
    # demography
    "add_landscape_change",
    "spDemography",
    # simulation
    "sim_ancestry",
    "sim_mutations",
    # plot
    "plot_model",
    "plot_landscape",
    "plot_timeseries",
    # analysis
    "filter_gt",
    "calc_sumstats",
]


# Mapping of public names to their source modules.  Submodules are imported
# lazily on first attribute access so that merely importing ``spaceprime``
# (e.g. when the CLI entry-point loads) does not pull in heavy scientific
# dependencies before they are actually needed.
_lazy_imports: dict[str, str] = {
    # utilities
    "create_raster": "spaceprime.utilities",
    "raster_to_demes": "spaceprime.utilities",
    "calc_migration_matrix": "spaceprime.utilities",
    "split_landscape_by_pop": "spaceprime.utilities",
    "mtp_thresh_from_coords": "spaceprime.utilities",
    "coords_to_sample_dict": "spaceprime.utilities",
    "anc_to_deme_dict": "spaceprime.utilities",
    "coords_to_deme_dict": "spaceprime.utilities",
    # demography
    "add_landscape_change": "spaceprime.demography",
    "spDemography": "spaceprime.demography",
    # simulation
    "sim_ancestry": "spaceprime.simulation",
    "sim_mutations": "spaceprime.simulation",
    # plot
    "plot_model": "spaceprime.plot",
    "plot_landscape": "spaceprime.plot",
    "plot_timeseries": "spaceprime.plot",
    # analysis
    "filter_gt": "spaceprime.analysis",
    "calc_sumstats": "spaceprime.analysis",
}


def __getattr__(name: str):
    if name in _lazy_imports:
        import importlib

        module = importlib.import_module(_lazy_imports[name])
        obj = getattr(module, name)
        # Cache in the module namespace so subsequent accesses are instant.
        globals()[name] = obj
        return obj
    raise AttributeError(f"module 'spaceprime' has no attribute {name!r}")


def __dir__():
    return list(globals()) + list(_lazy_imports)
