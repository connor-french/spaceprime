#!/usr/bin/env python

"""Tests for `cli` module."""
import numpy as np
import pytest

from spaceprime import demography, simulation
from spaceprime.cli import get_coal_times, get_map_dict, sci_notation_int


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def simple_demo():
    """2x2 spDemography with a single ancestral population."""
    d = np.array([[1000.0, 500.0], [200.0, 300.0]])
    demo = demography.spDemography()
    demo.stepping_stone_2d(d, rate=0.01)
    demo.add_ancestral_populations(anc_sizes=[5000], merge_time=500)
    return demo


@pytest.fixture
def raster_2x2(tmp_path):
    """Minimal 2x2 rasterio dataset matching the simple_demo grid."""
    import rasterio
    from rasterio.transform import from_bounds

    path = tmp_path / "raster.tif"
    transform = from_bounds(0, 0, 1, 1, 2, 2)
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=2,
        width=2,
        count=1,
        dtype="float32",
        crs="EPSG:4326",
        transform=transform,
    ) as dst:
        dst.write(np.ones((1, 2, 2), dtype="float32"))
    return rasterio.open(path)


# ---------------------------------------------------------------------------
# get_map_dict
# ---------------------------------------------------------------------------


class TestGetMapDict:
    def test_default_sample_num(self, simple_demo):
        """Default sample_num=2 is reflected in every deme's value."""
        result = get_map_dict(simple_demo)
        assert all(v == 2 for v in result.values())

    def test_custom_sample_num(self, simple_demo):
        """sample_num is passed through to each deme's value."""
        result = get_map_dict(simple_demo, sample_num=5)
        assert all(v == 5 for v in result.values())

    def test_excludes_small_demes(self, simple_demo):
        """Demes smaller than min_num_inds are excluded."""
        result_strict = get_map_dict(simple_demo, min_num_inds=400)
        result_loose = get_map_dict(simple_demo, min_num_inds=2)
        assert len(result_strict) < len(result_loose)


# ---------------------------------------------------------------------------
# get_coal_times
# ---------------------------------------------------------------------------


class TestGetCoalTimes:
    def test_output_shape_matches_raster(self, simple_demo, raster_2x2):
        """get_coal_times returns an array shaped like the raster."""
        ts = simulation.sim_ancestry(
            samples=get_map_dict(simple_demo, sample_num=2),
            demography=simple_demo,
            sequence_length=1000,
            random_seed=42,
        )
        ts = simulation.sim_mutations(ts, rate=1e-8, random_seed=42)
        result = get_coal_times(ts, raster_2x2, num_anc_pops=1, sample_num=2)
        assert result.shape == (2, 2)

    def test_custom_sample_num_changes_which_pops_are_scored(self, simple_demo, raster_2x2):
        """With sample_num=5, populations sampled at 2 get -1 (no match)."""
        ts_2 = simulation.sim_ancestry(
            samples=get_map_dict(simple_demo, sample_num=2),
            demography=simple_demo,
            sequence_length=1000,
            random_seed=42,
        )
        ts_2 = simulation.sim_mutations(ts_2, rate=1e-8, random_seed=42)
        result = get_coal_times(ts_2, raster_2x2, num_anc_pops=1, sample_num=5)
        # All demes were sampled at 2, not 5, so none should be scored
        assert np.all(result == -1)


def test_sci_notation_int_plain_int():
    assert sci_notation_int("1000") == 1000
    assert isinstance(sci_notation_int("1000"), int)


def test_sci_notation_int_scientific_notation():
    assert sci_notation_int("1e6") == 1000000
    assert isinstance(sci_notation_int("1e6"), int)


def test_sci_notation_int_scientific_notation_small_exponent():
    assert sci_notation_int("1e3") == 1000
    assert sci_notation_int("2.5e4") == 25000


def test_sci_notation_int_raises_on_non_numeric():
    with pytest.raises((ValueError, TypeError)):
        sci_notation_int("abc")
