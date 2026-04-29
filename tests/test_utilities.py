"""Tests for the `utilities` module."""

import numpy as np
import pytest

from spaceprime import utilities


# ---------------------------------------------------------------------------
# calc_migration_matrix
# ---------------------------------------------------------------------------


class TestCalcMigrationMatrix:
    """Tests for utilities.calc_migration_matrix."""

    def test_output_shape(self, deme_array_2d):
        """Migration matrix has shape (N, N) where N = rows * cols."""
        M = utilities.calc_migration_matrix(deme_array_2d, 0.01)
        n_demes = deme_array_2d.shape[0] * deme_array_2d.shape[1]
        assert M.shape == (n_demes, n_demes)

    def test_larger_grid_shape(self):
        """Shape scales correctly for a 3x4 deme grid."""
        demes = np.ones((3, 4)) * 500.0
        M = utilities.calc_migration_matrix(demes, 0.01, scale=False)
        assert M.shape == (12, 12)

    def test_diagonal_is_zero(self, deme_array_2d):
        """Self-migration (diagonal) is always zero."""
        M = utilities.calc_migration_matrix(deme_array_2d, 0.01)
        assert np.all(np.diag(M) == 0)

    def test_non_neighbors_are_zero(self, deme_array_2d):
        """Diagonal (non-adjacent) deme pairs have zero migration."""
        M = utilities.calc_migration_matrix(deme_array_2d, 0.01)
        # In a 2x2 grid deme_0_0 (index 0) and deme_1_1 (index 3) are diagonal
        assert M[0, 3] == 0
        assert M[3, 0] == 0

    def test_scaled_migration_asymmetric(self, deme_array_2d):
        """Scaled migration is donor/recipient * rate, so it is asymmetric."""
        M = utilities.calc_migration_matrix(deme_array_2d, 0.01, scale=True)
        # M[i, j] = (demes[j] / demes[i]) * rate
        # deme_0_0 (1000) -> deme_0_1 (500): (500/1000)*0.01 = 0.005
        assert np.isclose(M[0, 1], 0.005)
        # deme_0_1 (500) -> deme_0_0 (1000): (1000/500)*0.01 = 0.02
        assert np.isclose(M[1, 0], 0.02)

    def test_constant_migration_symmetric(self, deme_array_2d):
        """With scale=False all adjacent non-zero deme pairs get the same rate."""
        M = utilities.calc_migration_matrix(deme_array_2d, 0.01, scale=False)
        assert np.isclose(M[0, 1], 0.01)
        assert np.isclose(M[1, 0], 0.01)
        assert np.isclose(M[0, 2], 0.01)
        assert np.isclose(M[2, 0], 0.01)

    def test_zero_population_deme_gets_zero_migration(self):
        """Any migration involving an empty deme is set to zero."""
        demes = np.array([[1000.0, 1e-10], [200.0, 300.0]])
        M = utilities.calc_migration_matrix(demes, 0.01)
        # deme_0_1 (index 1) is effectively empty
        assert M[0, 1] == 0
        assert M[1, 0] == 0
        assert M[1, 3] == 0
        assert M[3, 1] == 0

    def test_type_error_demes_not_array(self):
        """Passing a list instead of ndarray raises TypeError."""
        with pytest.raises(TypeError, match="numpy array"):
            utilities.calc_migration_matrix([[1.0, 2.0], [3.0, 4.0]], 0.01)

    def test_type_error_rate_not_float(self):
        """Passing an int rate raises TypeError."""
        with pytest.raises(TypeError, match="float"):
            utilities.calc_migration_matrix(np.array([[1000.0, 500.0]]), 1)

    def test_type_error_scale_not_bool(self):
        """Passing a string for scale raises TypeError."""
        with pytest.raises(TypeError, match="bool"):
            utilities.calc_migration_matrix(
                np.array([[1000.0, 500.0]]), 0.01, scale="yes"
            )


# ---------------------------------------------------------------------------
# raster_to_demes
# ---------------------------------------------------------------------------


class TestRasterToDemes:
    """Tests for utilities.raster_to_demes."""

    @pytest.mark.parametrize("transformation", ["linear", "threshold", "sigmoid"])
    def test_returns_ndarray(self, transformation):
        """Each transformation mode returns a numpy ndarray."""
        arr = np.array([[0.8, 0.4], [0.2, 0.6]])
        kwargs = {"threshold": 0.5} if transformation == "threshold" else {}
        result = utilities.raster_to_demes(
            arr, transformation=transformation, max_local_size=1000, **kwargs
        )
        assert isinstance(result, np.ndarray)

    def test_linear_scales_by_max_local_size(self):
        """Linear transformation multiplies values by max_local_size."""
        arr = np.array([[1.0, 0.5], [0.25, 0.0]])
        result = utilities.raster_to_demes(
            arr, transformation="linear", max_local_size=1000
        )
        assert result[0, 0] == 1000
        assert result[0, 1] == 500
        assert result[1, 1] <= 1e-9  # 0.0 maps to effectively zero

    def test_linear_threshold_zeroes_low_values(self):
        """Linear mode with threshold=0.5 zeroes cells below 0.5*max_local_size."""
        arr = np.array([[0.8, 0.3], [0.1, 0.9]])
        result = utilities.raster_to_demes(
            arr, transformation="linear", max_local_size=1000, threshold=0.5
        )
        assert result[0, 0] > 1e-9   # 0.8 >= threshold
        assert result[1, 1] > 1e-9   # 0.9 >= threshold
        assert result[0, 1] <= 1e-9  # 0.3 < threshold
        assert result[1, 0] <= 1e-9  # 0.1 < threshold

    def test_threshold_transformation_binary(self):
        """Threshold mode produces only max_local_size or ~zero values."""
        arr = np.array([[0.8, 0.4], [0.2, 0.6]])
        result = utilities.raster_to_demes(
            arr, transformation="threshold", max_local_size=1000, threshold=0.5
        )
        assert result[0, 0] == 1000   # 0.8 >= 0.5
        assert result[1, 1] == 1000   # 0.6 >= 0.5
        assert result[0, 1] <= 1e-9  # 0.4 < 0.5
        assert result[1, 0] <= 1e-9  # 0.2 < 0.5

    def test_sigmoid_transformation_bounded(self):
        """Sigmoid output is bounded by max_local_size and all positive."""
        arr = np.array([[0.8, 0.4], [0.2, 0.6]])
        result = utilities.raster_to_demes(
            arr,
            transformation="sigmoid",
            max_local_size=1000,
            inflection_point=0.5,
            slope=0.1,
        )
        assert np.all(result > 0)
        # Values above inflection point should be larger than those below
        assert result[0, 0] > result[0, 1]  # 0.8 > 0.4

    def test_normalize_flag_scales_to_0_1(self):
        """normalize=True rescales input so max maps to max_local_size."""
        arr = np.array([[100.0, 200.0], [300.0, 400.0]])
        result = utilities.raster_to_demes(
            arr, transformation="linear", max_local_size=1000, normalize=True
        )
        assert result.max() == 1000

    def test_rasterio_input_accepted(self, raster):
        """A rasterio DatasetReader object is accepted as input."""
        result = utilities.raster_to_demes(
            raster, transformation="linear", max_local_size=1000
        )
        assert isinstance(result, np.ndarray)
        assert result.ndim == 3  # (bands, rows, cols)

    def test_type_error_invalid_raster(self):
        """Passing a string raises TypeError."""
        with pytest.raises(TypeError, match="numpy array or rasterio"):
            utilities.raster_to_demes("not_a_raster", transformation="linear")

    def test_value_error_invalid_transformation(self):
        """Passing an unknown transformation name raises ValueError."""
        arr = np.array([[0.5, 0.3]])
        with pytest.raises(ValueError, match="Invalid transformation"):
            utilities.raster_to_demes(arr, transformation="invalid")


# ---------------------------------------------------------------------------
# anc_to_deme_dict
# ---------------------------------------------------------------------------


class TestAncToDemeDict:
    """Tests for utilities.anc_to_deme_dict."""

    def test_basic_mapping(self):
        """Deme indices are grouped correctly by ancestral population ID."""
        anc_pops = np.array([[1, 1], [2, 2]])
        deme_dict = {0: 2, 1: 4, 2: 3, 3: 2}
        result = utilities.anc_to_deme_dict(anc_pops, deme_dict)
        assert sorted(result[1]) == [0, 1]
        assert sorted(result[2]) == [2, 3]

    def test_one_deme_per_ancestral_pop(self):
        """Works when every deme maps to a distinct ancestral population."""
        anc_pops = np.array([[1, 2], [3, 4]])
        deme_dict = {0: 2, 1: 2, 2: 2, 3: 2}
        result = utilities.anc_to_deme_dict(anc_pops, deme_dict)
        assert len(result) == 4
        for key in [1, 2, 3, 4]:
            assert len(result[key]) == 1
