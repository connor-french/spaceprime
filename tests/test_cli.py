#!/usr/bin/env python

"""Tests for `cli` module."""
import argparse

import numpy as np
import pytest
from geopandas import GeoDataFrame
from shapely.geometry import Point

from spaceprime import demography, simulation
from spaceprime.cli import (
    build_parser,
    get_coal_times,
    get_map_dict,
    read_individuals,
    sci_notation_int,
    setup_demography,
)


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


# ---------------------------------------------------------------------------
# read_individuals
# ---------------------------------------------------------------------------


class TestReadIndividuals:
    def test_returns_none_when_individuals_omitted(self):
        """No UnboundLocalError when --individuals is not supplied."""
        args = argparse.Namespace(individuals=None)
        assert read_individuals(args, coords=[]) is None

    def test_returns_list_when_given_list(self):
        """Inline list of IDs is passed through unchanged."""
        ids = ["ind_1", "ind_2"]
        args = argparse.Namespace(individuals=ids)
        coords = [None, None]
        assert read_individuals(args, coords) == ids

    def test_raises_when_length_mismatch(self):
        """Mismatched individuals / coords length raises ValueError."""
        args = argparse.Namespace(individuals=["ind_1", "ind_2", "ind_3"])
        coords = [None, None]
        with pytest.raises(ValueError, match="Number of individuals"):
            read_individuals(args, coords)


# ---------------------------------------------------------------------------
# Argument parser defaults and type lambdas
# ---------------------------------------------------------------------------


class TestBuildParser:
    def test_mig_rate_default_is_1e8(self):
        """--mig_rate defaults to [1e-8], not None."""
        args = build_parser().parse_args([])
        assert args.mig_rate == [1e-8]

    def test_anc_sizes_single_value_is_scalar_int(self):
        """Single --anc_sizes token parses as a plain int, not a list."""
        args = build_parser().parse_args(["-as", "5000"])
        assert args.anc_sizes == [5000]
        assert isinstance(args.anc_sizes[0], int)

    def test_anc_sizes_comma_pair_is_list(self):
        """Comma-separated --anc_sizes pair parses as a list of two ints."""
        args = build_parser().parse_args(["-as", "5000,8000"])
        assert args.anc_sizes == [[5000, 8000]]

    def test_anc_sizes_multiple_singles(self):
        """Multiple space-separated --anc_sizes values parse as a flat int list."""
        args = build_parser().parse_args(["-as", "5000", "6000"])
        assert args.anc_sizes == [5000, 6000]
        assert all(isinstance(v, int) for v in args.anc_sizes)

    def test_anc_sizes_defaults_to_none(self):
        """--anc_sizes defaults to None when omitted."""
        args = build_parser().parse_args([])
        assert args.anc_sizes is None


# ---------------------------------------------------------------------------
# setup_demography
# ---------------------------------------------------------------------------


class TestSetupDemography:
    def test_returns_spdemography_with_default_mig_rate(self, raster):
        """setup_demography succeeds with mig_rate=0.0 and no ancestral pops."""
        result = setup_demography(
            raster=raster,
            coords=None,
            max_local_size=1000,
            threshold=None,
            inflection_point=0.5,
            slope=0.05,
            mig_rate=0.0,
            scale=True,
            anc_pop_id=None,
            timesteps=1,
            anc_sizes=None,
            merge_time=None,
            anc_merge_time=None,
            anc_merge_size=None,
            anc_mig_rate=None,
        )
        assert isinstance(result, demography.spDemography)
        assert len(result.populations) > 0

    def test_anc_sizes_none_skips_ancestral_populations(self, raster):
        """No ancestral populations are added when anc_sizes is None."""
        result = setup_demography(
            raster=raster,
            coords=None,
            max_local_size=1000,
            threshold=None,
            inflection_point=0.5,
            slope=0.05,
            mig_rate=0.0,
            scale=True,
            anc_pop_id=None,
            timesteps=1,
            anc_sizes=None,
            merge_time=None,
            anc_merge_time=None,
            anc_merge_size=None,
            anc_mig_rate=None,
        )
        assert not any("ANC" in p.name for p in result.populations)

    def test_anc_sizes_provided_adds_ancestral_population(self, raster):
        """Ancestral population is added when anc_sizes and merge_time are set."""
        result = setup_demography(
            raster=raster,
            coords=None,
            max_local_size=1000,
            threshold=None,
            inflection_point=0.5,
            slope=0.05,
            mig_rate=0.0,
            scale=True,
            anc_pop_id=None,
            timesteps=1,
            anc_sizes=[5000],
            merge_time=500,
            anc_merge_time=None,
            anc_merge_size=None,
            anc_mig_rate=None,
        )
        assert any("ANC" in p.name for p in result.populations)
