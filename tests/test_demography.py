"""Tests for the `demography` module."""

import numpy as np
import pytest

from spaceprime import demography


# ---------------------------------------------------------------------------
# spDemography.stepping_stone_2d
# ---------------------------------------------------------------------------


class TestSteppingStone2D:
    """Tests for spDemography.stepping_stone_2d."""

    def test_population_count(self, deme_array_2d):
        """Number of populations equals the total number of deme cells."""
        demo = demography.spDemography()
        demo.stepping_stone_2d(deme_array_2d, rate=0.01)
        n_demes = deme_array_2d.shape[0] * deme_array_2d.shape[1]
        assert len(demo.populations) == n_demes

    def test_population_names_follow_convention(self, deme_array_2d):
        """Populations are named deme_i_j for each row i, column j."""
        demo = demography.spDemography()
        demo.stepping_stone_2d(deme_array_2d, rate=0.01)
        names = [p.name for p in demo.populations]
        assert "deme_0_0" in names
        assert "deme_0_1" in names
        assert "deme_1_0" in names
        assert "deme_1_1" in names

    def test_demes_attribute_set(self, deme_array_2d):
        """The `demes` attribute is set to the input array."""
        demo = demography.spDemography()
        demo.stepping_stone_2d(deme_array_2d, rate=0.01)
        assert np.array_equal(demo.demes, deme_array_2d)

    def test_migration_array_shape_2d_input(self, deme_array_2d):
        """With a 2D deme array, migration_array is a 2D (N×N) matrix."""
        demo = demography.spDemography()
        demo.stepping_stone_2d(deme_array_2d, rate=0.01)
        n_demes = deme_array_2d.shape[0] * deme_array_2d.shape[1]
        assert demo.migration_array.shape == (n_demes, n_demes)

    def test_scale_true_produces_asymmetric_rates(self, deme_array_2d):
        """Scaled migration yields different rates in each direction."""
        demo = demography.spDemography()
        demo.stepping_stone_2d(deme_array_2d, rate=0.01, scale=True)
        M = demo.migration_array
        # deme_0_0 <-> deme_0_1: sizes 1000 vs 500, so rates differ
        assert not np.isclose(M[0, 1], M[1, 0])

    def test_scale_false_produces_constant_rates(self, deme_array_2d):
        """With scale=False all adjacent populated deme pairs share the same rate."""
        demo = demography.spDemography()
        demo.stepping_stone_2d(deme_array_2d, rate=0.01, scale=False)
        M = demo.migration_array
        assert np.isclose(M[0, 1], 0.01)
        assert np.isclose(M[1, 0], 0.01)

    def test_3d_input_migration_array_is_3d(self, deme_array_3d):
        """With a 3D (T×R×C) deme array, migration_array has shape (T, N, N)."""
        demo = demography.spDemography()
        demo.stepping_stone_2d(deme_array_3d, rate=0.01, timesteps=100)
        T = deme_array_3d.shape[0]
        N = deme_array_3d.shape[1] * deme_array_3d.shape[2]
        assert demo.migration_array.shape == (T, N, N)

    def test_3d_input_demes_attribute_shape(self, deme_array_3d):
        """The `demes` attribute preserves the 3D shape of the input."""
        demo = demography.spDemography()
        demo.stepping_stone_2d(deme_array_3d, rate=0.01, timesteps=100)
        assert demo.demes.shape == deme_array_3d.shape

    def test_custom_2d_migration_matrix_accepted(self, deme_array_2d):
        """A custom N×N migration matrix is stored without modification."""
        n = deme_array_2d.shape[0] * deme_array_2d.shape[1]
        custom_mig = np.full((n, n), 0.005)
        np.fill_diagonal(custom_mig, 0.0)
        demo = demography.spDemography()
        demo.stepping_stone_2d(deme_array_2d, rate=custom_mig)
        assert np.allclose(demo.migration_array, custom_mig)

    def test_custom_migration_wrong_shape_raises(self, deme_array_2d):
        """A custom migration matrix with wrong shape raises AssertionError."""
        wrong_shape = np.ones((3, 3))  # 2x2 grid needs 4x4
        demo = demography.spDemography()
        with pytest.raises(AssertionError):
            demo.stepping_stone_2d(deme_array_2d, rate=wrong_shape)


# ---------------------------------------------------------------------------
# spDemography.add_ancestral_populations
# ---------------------------------------------------------------------------


class TestAddAncestralPopulations:
    """Tests for spDemography.add_ancestral_populations."""

    def test_single_ancestral_population_added(self, basic_demo):
        """One ANC population is added when anc_id is not provided."""
        basic_demo.add_ancestral_populations(anc_sizes=[5000], merge_time=500)
        anc_pops = [p for p in basic_demo.populations if "ANC" in p.name]
        assert len(anc_pops) == 1
        assert anc_pops[0].name == "ANC_1"

    def test_single_ancestral_population_size(self, basic_demo):
        """The ancestral population is created with the specified initial size."""
        basic_demo.add_ancestral_populations(anc_sizes=[5000], merge_time=500)
        anc_pop = next(p for p in basic_demo.populations if p.name == "ANC_1")
        assert anc_pop.initial_size == 5000

    def test_multiple_ancestral_populations(self, basic_demo):
        """Two ANC populations are added when anc_id assigns demes to two groups."""
        anc_id = np.array([[1, 1], [2, 2]])
        basic_demo.add_ancestral_populations(
            anc_sizes=[3000, 2000], merge_time=500, anc_id=anc_id
        )
        anc_pops = [p for p in basic_demo.populations if "ANC" in p.name]
        anc_names = {p.name for p in anc_pops}
        assert "ANC_1" in anc_names
        assert "ANC_2" in anc_names

    def test_duplicate_add_raises_value_error(self, basic_demo):
        """Adding ancestral populations twice raises ValueError."""
        basic_demo.add_ancestral_populations(anc_sizes=[5000], merge_time=500)
        with pytest.raises(ValueError, match="already contains ancestral"):
            basic_demo.add_ancestral_populations(anc_sizes=[5000], merge_time=500)

    def test_list_merge_time_raises_value_error(self, basic_demo):
        """Passing a list as merge_time raises ValueError."""
        with pytest.raises(ValueError, match="single float or int"):
            basic_demo.add_ancestral_populations(
                anc_sizes=[5000], merge_time=[100, 200]
            )

    def test_non_numeric_merge_time_raises_type_error(self, basic_demo):
        """Passing a string as merge_time raises TypeError."""
        with pytest.raises(TypeError, match="float or int"):
            basic_demo.add_ancestral_populations(
                anc_sizes=[5000], merge_time="large"
            )

    def test_migration_rate_between_ancestral_pops(self, basic_demo):
        """When migration_rate is given, ANC populations have nonzero migration."""
        anc_id = np.array([[1, 1], [2, 2]])
        basic_demo.add_ancestral_populations(
            anc_sizes=[3000, 2000],
            merge_time=500,
            anc_id=anc_id,
            migration_rate=0.001,
        )
        # The migration matrix should have grown to include ANC populations
        n_total = len(basic_demo.populations)
        assert len(basic_demo.migration_matrix) == n_total

    def test_scalar_anc_sizes_converted_to_list(self, basic_demo):
        """Passing a scalar for anc_sizes does not raise an error."""
        basic_demo.add_ancestral_populations(anc_sizes=5000, merge_time=500)
        anc_pops = [p for p in basic_demo.populations if "ANC" in p.name]
        assert len(anc_pops) == 1
