"""Smoke tests for the top-level `spaceprime` package."""

import spaceprime


def test_package_exposes_spdemography():
    assert hasattr(spaceprime, "spDemography")


def test_package_exposes_simulation_wrappers():
    assert hasattr(spaceprime, "sim_ancestry")
    assert hasattr(spaceprime, "sim_mutations")


def test_package_exposes_utility_functions():
    for name in (
        "raster_to_demes",
        "calc_migration_matrix",
        "coords_to_sample_dict",
        "anc_to_deme_dict",
        "coords_to_deme_dict",
        "split_landscape_by_pop",
        "mtp_thresh_from_coords",
        "create_raster",
    ):
        assert hasattr(spaceprime, name), f"spaceprime.{name} not found"


def test_package_exposes_plot_functions():
    for name in ("plot_model", "plot_landscape", "plot_timeseries"):
        assert hasattr(spaceprime, name), f"spaceprime.{name} not found"


def test_package_exposes_analysis_functions():
    assert hasattr(spaceprime, "filter_gt")
    assert hasattr(spaceprime, "calc_sumstats")
    # The analysis extras may be unavailable in minimal installs
    if spaceprime.filter_gt is not None:
        assert callable(spaceprime.filter_gt)
