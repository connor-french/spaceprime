"""Shared pytest fixtures for the spaceprime test suite."""

import numpy as np
import pytest
import rasterio

from spaceprime import demography


@pytest.fixture
def deme_array_2d():
    """A simple 2x2 array of deme sizes."""
    return np.array([[1000.0, 500.0], [200.0, 300.0]])


@pytest.fixture
def deme_array_3d():
    """A 2-timestep, 2x2 array of deme sizes (most-recent first)."""
    return np.array(
        [
            [[1000.0, 500.0], [200.0, 300.0]],
            [[800.0, 400.0], [100.0, 200.0]],
        ]
    )


@pytest.fixture
def raster():
    """A 2x2 rasterio DatasetReader (EPSG:4326, cell centres at 0.5/-0.5 etc.)."""
    return rasterio.open("tests/data/demes_raster.tif")


@pytest.fixture
def basic_demo(deme_array_2d):
    """An spDemography built from the 2x2 deme array with scaled migration."""
    demo = demography.spDemography()
    demo.stepping_stone_2d(deme_array_2d, rate=0.01)
    return demo


@pytest.fixture
def demo_with_anc(deme_array_2d):
    """An spDemography with a single ancestral population added."""
    demo = demography.spDemography()
    demo.stepping_stone_2d(deme_array_2d, rate=0.01)
    demo.add_ancestral_populations(anc_sizes=[5000], merge_time=500)
    return demo
