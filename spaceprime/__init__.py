"""
spaceprime is a python package to facilitate spatially explicit coalescent modeling in msprime.
Please see the documentation at https://connor-french.github.io/spaceprime
for more information.
"""

__author__ = """Connor French"""
__email__ = "french.connor.m@gmail.com"
__version__ = "0.0.4"


from spaceprime.utilities import (
    create_raster,
    raster_to_demes,
    calc_migration_matrix,
    split_landscape_by_pop,
    mtp_thresh_from_coords,
    coords_to_sample_dict,
    anc_to_deme_dict,
    coords_to_deme_dict,
)

from spaceprime.demography import (
    add_landscape_change,
    spDemography,
)

from spaceprime.simulation import (
    sim_ancestry,
    sim_mutations,
)

from spaceprime.plot import (
    plot_model,
    plot_landscape,
    plot_timeseries,
)

# optional install of analysis tools
try:
    from spaceprime.analysis import filter_gt, calc_sumstats
except ImportError:
    filter_gt = None
    calc_sumstats = None
