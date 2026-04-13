"""Tests for the `analysis` module."""

import numpy as np
import pandas as pd
import pytest

from spaceprime import analysis


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def random_gt():
    """A random haploid genotype matrix: 200 loci × 32 haplotypes (16 diploids)."""
    np.random.seed(42)
    return np.random.randint(0, 2, size=(200, 32))


@pytest.fixture
def deme_dict_inds():
    """Indices of 4 individuals in each of 4 demes (diploid, 0-indexed)."""
    return {k: list(range(k * 4, (k + 1) * 4)) for k in range(4)}


@pytest.fixture
def grid_deme_dict_inds():
    """10 demes × 4 individuals — used for spatial statistics (needs >4 demes)."""
    return {k: list(range(k * 4, (k + 1) * 4)) for k in range(10)}


@pytest.fixture
def grid_coords_dict():
    """2×5 grid of lon/lat coordinates (non-collinear, needed for Voronoi)."""
    return {
        0: [0.5, -0.5],
        1: [1.5, -0.5],
        2: [2.5, -0.5],
        3: [3.5, -0.5],
        4: [4.5, -0.5],
        5: [0.5, -2.5],
        6: [1.5, -2.5],
        7: [2.5, -2.5],
        8: [3.5, -2.5],
        9: [4.5, -2.5],
    }


@pytest.fixture
def filtered_outputs(deme_dict_inds, random_gt):
    """Pre-filtered genotype data for a 4-deme case."""
    np.random.seed(0)
    return analysis.filter_gt(random_gt, deme_dict_inds=deme_dict_inds)


@pytest.fixture
def grid_filtered_outputs(grid_deme_dict_inds):
    """Pre-filtered genotype data for the 10-deme spatial statistics case."""
    np.random.seed(0)
    gt = np.random.randint(0, 2, size=(500, 80))
    return analysis.filter_gt(gt, deme_dict_inds=grid_deme_dict_inds)


# ---------------------------------------------------------------------------
# filter_gt
# ---------------------------------------------------------------------------


class TestFilterGt:
    """Tests for analysis.filter_gt."""

    def test_returns_four_element_tuple(self, random_gt):
        """filter_gt always returns a 4-tuple."""
        result = analysis.filter_gt(random_gt)
        assert len(result) == 4

    def test_output_genotype_array_type(self, random_gt):
        """First output is a scikit-allel GenotypeArray."""
        import allel

        gt_out, _, _, _ = analysis.filter_gt(random_gt)
        assert isinstance(gt_out, allel.GenotypeArray)

    def test_output_allele_counts_type(self, random_gt):
        """Second output is a scikit-allel AlleleCountsArray."""
        import allel

        _, ac_out, _, _ = analysis.filter_gt(random_gt)
        assert isinstance(ac_out, allel.AlleleCountsArray)

    def test_no_deme_dict_returns_none_ac_demes(self, random_gt):
        """ac_demes is None when deme_dict_inds is not provided."""
        _, _, ac_demes, _ = analysis.filter_gt(random_gt)
        assert ac_demes is None

    def test_deme_dict_returns_per_deme_allele_counts(
        self, random_gt, deme_dict_inds
    ):
        """When deme_dict_inds is given, ac_demes contains one entry per deme."""
        _, _, ac_demes, _ = analysis.filter_gt(
            random_gt, deme_dict_inds=deme_dict_inds
        )
        assert ac_demes is not None
        assert set(ac_demes.keys()) == set(deme_dict_inds.keys())

    def test_filter_monomorphic_reduces_loci(self, random_gt):
        """filter_monomorphic=True produces fewer loci than filter_monomorphic=False."""
        np.random.seed(1)
        gt = np.random.randint(0, 2, size=(200, 32))
        # Force the first 10 rows to be monomorphic
        gt[:10, :] = 0

        gt_filtered, _, _, _ = analysis.filter_gt(
            gt, filter_monomorphic=True, filter_singletons=False
        )
        gt_unfiltered, _, _, _ = analysis.filter_gt(
            gt, filter_monomorphic=False, filter_singletons=False
        )
        assert gt_filtered.shape[0] <= gt_unfiltered.shape[0]

    def test_filter_singletons_removes_singletons(self, random_gt):
        """After filtering, no singleton sites remain in the allele counts."""
        import allel

        np.random.seed(2)
        _, ac_out, _, _ = analysis.filter_gt(random_gt, filter_singletons=True)
        if ac_out.shape[0] > 0:
            assert not ac_out.is_singleton(allele=1).any()

    def test_filter_singletons_false_retains_more_sites(self, random_gt):
        """Disabling singleton filtering keeps at least as many loci."""
        np.random.seed(3)
        gt_no_single, _, _, _ = analysis.filter_gt(
            random_gt, filter_singletons=True
        )
        gt_with_single, _, _, _ = analysis.filter_gt(
            random_gt, filter_singletons=False
        )
        assert gt_no_single.shape[0] <= gt_with_single.shape[0]

    def test_missing_data_perc_does_not_crash(self, random_gt, deme_dict_inds):
        """filter_gt runs without error when missing_data_perc > 0."""
        np.random.seed(4)
        result = analysis.filter_gt(
            random_gt,
            deme_dict_inds=deme_dict_inds,
            missing_data_perc=0.2,
        )
        assert len(result) == 4


# ---------------------------------------------------------------------------
# calc_sumstats
# ---------------------------------------------------------------------------


class TestCalcSumstats:
    """Tests for analysis.calc_sumstats."""

    def test_returns_dict_without_ac_demes(self, filtered_outputs):
        """Without ac_demes, calc_sumstats returns a plain dict."""
        gt_out, ac_out, _, _ = filtered_outputs
        coords_dict = {0: [0.5, -0.5]}
        result = analysis.calc_sumstats(ac_out, coords_dict)
        assert isinstance(result, dict)

    def test_global_summary_stat_keys_present(self, filtered_outputs):
        """Standard global statistics keys are present in the output."""
        _, ac_out, _, _ = filtered_outputs
        coords_dict = {0: [0.5, -0.5]}
        result = analysis.calc_sumstats(ac_out, coords_dict)
        for key in ("pi", "pi_sd", "sfs_h1", "sfs_h2"):
            assert key in result, f"Expected key '{key}' missing from output"

    def test_return_df_true_gives_dataframe(self, filtered_outputs):
        """return_df=True returns a pandas DataFrame with one row."""
        _, ac_out, _, _ = filtered_outputs
        coords_dict = {0: [0.5, -0.5]}
        result = analysis.calc_sumstats(ac_out, coords_dict, return_df=True)
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 1

    def test_precision_rounds_values(self, filtered_outputs):
        """Non-NaN values are rounded to the specified number of decimal places."""
        _, ac_out, _, _ = filtered_outputs
        coords_dict = {0: [0.5, -0.5]}
        result = analysis.calc_sumstats(ac_out, coords_dict, precision=3)
        for val in result.values():
            f = float(val)
            if np.isnan(f):
                continue  # NaN is valid (e.g. Tajima's D with too few loci)
            assert np.isclose(f, round(f, 3), atol=1e-9)

    def test_per_deme_pi_keys_added_with_ac_demes(
        self, grid_filtered_outputs, grid_coords_dict, grid_deme_dict_inds
    ):
        """pi_deme_N keys are present for each deme when ac_demes is given."""
        _, ac_out, ac_demes, _ = grid_filtered_outputs
        result = analysis.calc_sumstats(
            ac_out, grid_coords_dict, ac_demes=ac_demes
        )
        for deme_id in grid_deme_dict_inds:
            assert f"pi_deme_{deme_id}" in result

    def test_spatial_stats_added_with_ac_demes(
        self, grid_filtered_outputs, grid_coords_dict
    ):
        """Moran's I and IBD slope/R² are included when ac_demes is supplied."""
        _, ac_out, ac_demes, _ = grid_filtered_outputs
        result = analysis.calc_sumstats(
            ac_out, grid_coords_dict, ac_demes=ac_demes
        )
        assert "morans_i" in result
        assert "ibd_slope" in result
        assert "ibd_r2" in result

    def test_dxy_keys_added_for_all_deme_pairs(
        self, grid_filtered_outputs, grid_coords_dict, grid_deme_dict_inds
    ):
        """A dxy_i_j key is present for every pair of demes."""
        import itertools

        _, ac_out, ac_demes, _ = grid_filtered_outputs
        result = analysis.calc_sumstats(
            ac_out, grid_coords_dict, ac_demes=ac_demes
        )
        deme_ids = list(grid_deme_dict_inds.keys())
        for i, j in itertools.combinations(deme_ids, 2):
            assert f"dxy_{i}_{j}" in result
