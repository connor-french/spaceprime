"""
spaceprime is a python package to facilitate spatially explicit coalescent modeling in msprime.
Please see the documentation at https://connor-french.github.io/spaceprime
for more information.
"""

__author__ = """Connor French"""
__email__ = "french.connor.m@gmail.com"
__version__ = "0.0.2"


from spaceprime.utilities import (
    create_raster,
    raster_to_demes,
    calc_migration_matrix,
    split_landscape_by_pop,
    max_thresh_from_coords,
    coords_to_sample_dict,
    anc_to_deme_dict,
    sampled_cells_to_coords,
)
