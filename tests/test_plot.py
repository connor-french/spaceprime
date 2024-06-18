#!/usr/bin/env python

"""Tests for `plot` module."""
import pytest
import numpy as np
from spaceprime import plot
from spaceprime import demography
import rasterio
from math import isclose
import matplotlib


@pytest.fixture
def demes():
    return np.array([[1000, 500], [200, 300]])


@pytest.fixture
def mig_mat():
    # 0.01 is the migration rate between demes with scaling
    return np.array(
        [
            [0.0, 0.005, 0.002, 0.0],
            [0.02, 0.0, 0.0, 0.006],
            [0.05, 0.0, 0.0, 0.015],
            [0.0, 0.01666667, 0.00666667, 0.0],
        ]
    )


@pytest.fixture
def demo(demes, mig_mat):
    d = demography.spDemography()
    # Add demes to the demography object
    d.add_population(name="deme1", initial_size=demes[0][0])
    d.add_population(name="deme2", initial_size=demes[0][1])
    d.add_population(name="deme3", initial_size=demes[1][0])
    d.add_population(name="deme4", initial_size=demes[1][1])
    # Add migration rates between demes
    d.set_migration_rate("deme1", "deme2", rate=mig_mat[0][1])
    d.set_migration_rate("deme1", "deme3", rate=mig_mat[0][2])
    d.set_migration_rate("deme1", "deme4", rate=mig_mat[0][3])
    d.set_migration_rate("deme2", "deme1", rate=mig_mat[1][0])
    d.set_migration_rate("deme2", "deme3", rate=mig_mat[1][2])
    d.set_migration_rate("deme2", "deme4", rate=mig_mat[1][3])
    d.set_migration_rate("deme3", "deme1", rate=mig_mat[2][0])
    d.set_migration_rate("deme3", "deme2", rate=mig_mat[2][1])
    d.set_migration_rate("deme3", "deme4", rate=mig_mat[2][3])
    d.set_migration_rate("deme4", "deme1", rate=mig_mat[3][0])
    d.set_migration_rate("deme4", "deme2", rate=mig_mat[3][1])
    d.set_migration_rate("deme4", "deme3", rate=mig_mat[3][2])

    # add demes and migration_array attributes to the spDemography object
    d.demes = []
    d.demes.append(demes)

    d.migration_array = []
    d.migration_array.append(mig_mat)

    return d


@pytest.fixture
def raster():
    return rasterio.open("tests/data/demes_raster.tif")


# def write_rasterio_dataset(demes):
#   # Create a raster with the same values as demes
#   raster_shape = demes.shape
#   raster_transform = rasterio.transform.from_origin(0, 0, 1, 1)  # Assuming pixel size of 1
#   raster_crs = rasterio.crs.CRS.from_epsg(4326)  # Assuming EPSG code 4326 for geographic coordinates

#   with rasterio.open("tests/data/demes_raster.tif", "w", driver="GTiff", height=raster_shape[0], width=raster_shape[1],
#               count=1, dtype=demes.dtype, crs=raster_crs, transform=raster_transform) as dst:
#     dst.write(demes, 1)  # Write the raster values to the dataset


def test_get_outgoing_migration_rates_values(demes, mig_mat):
    outgoing_migration_rates = plot.get_outgoing_migration_rates(demes, mig_mat)
    assert outgoing_migration_rates.shape == (4, 5)
    assert outgoing_migration_rates["deme_size"].sum() == 2000
    assert isclose(
        outgoing_migration_rates["mig_north"].sum(), 0.066666667, rel_tol=1e-5
    )
    assert outgoing_migration_rates["mig_east"].sum() == 0.02
    assert outgoing_migration_rates["mig_south"].sum() == 0.008
    assert isclose(outgoing_migration_rates["mig_west"].sum(), 0.02666666, rel_tol=1e-5)


def test_plot_model_timestep(demo, raster) -> None:
    # check to make sure the function raises a ValueError if the timestep is invalid
    with pytest.raises(ValueError):
        plot.plot_model(demo, raster, 2)


def test_plot_model_demo(demo, raster) -> None:
    d1 = demography.spDemography()
    with pytest.raises(ValueError):
        plot.plot_model(d1, raster, 0)


def test_plot_landscape_demo(demo, raster):
    d1 = demography.spDemography()
    with pytest.raises(ValueError):
        plot.plot_landscape(d1, raster, 0)


def test_plot_landscape_timestep(demo, raster):
    with pytest.raises(ValueError):
        plot.plot_landscape(demo, raster, 2)


def test_plot_landscape_returns_plot(demo, raster):
    p = plot.plot_landscape(demo, raster, 0)
    assert isinstance(p, matplotlib.axes.Axes)
