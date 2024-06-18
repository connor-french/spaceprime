"""
spaceprime is a python package to facilitate spatially explicit coalescent modeling in msprime.
Please see the documentation at https://connor-french.github.io/spaceprime
for more information.
"""

__author__ = """Connor French"""
__email__ = "french.connor.m@gmail.com"
__version__ = "0.0.3"


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

# optional install of analysis tools
try:
    from spaceprime.analysis import filter_gt, calc_sumstats
except ImportError:
    filter_gt = None
    calc_sumstats = None

# optional install of plotting tools
try:
    from spaceprime.plot import get_outgoing_migration_rates, plot_model
except ImportError:
    get_outgoing_migration_rates = None
    plot_model = None
