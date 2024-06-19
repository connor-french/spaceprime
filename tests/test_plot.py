#!/usr/bin/env python

"""Tests for `plot` module."""
import pytest
import numpy as np
from spaceprime import plot
from spaceprime import demography
import rasterio
from math import isclose
import matplotlib
from matplotlib import pyplot as plt
import folium


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

    # add a second timestep to the demes array
    d.demes.append(demes * 2)

    d.migration_array = []
    d.migration_array.append(mig_mat)

    return d


@pytest.fixture
def raster():
    return rasterio.open("tests/data/demes_raster.tif")


# get_outgoing_migration_rates
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


def test_get_outgoing_migration_rates_navalue(demes, mig_mat):
    """Test that get_outgoing_migration_rates skips NA values in the migration matrix"""

    omr = plot.get_outgoing_migration_rates(demes, mig_mat, na_value=200)

    assert omr.shape == (3, 5)


# plot_model
def test_plot_model_timestep(demo, raster) -> None:
    """Test that the plot_model function raises a ValueError if the timestep is invalid"""

    with pytest.raises(ValueError):
        plot.plot_model(demo, raster, 2)


def test_plot_model_demo(demo, raster) -> None:
    """Test that the plot_model function raises a ValueError if the demography object isn't populated with deme sizes yet"""

    d1 = demography.spDemography()
    with pytest.raises(ValueError):
        plot.plot_model(d1, raster, 0)


def test_plot_model_returns_plot(demo, raster) -> None:
    # Call the function
    p = plot.plot_model(demo, raster, 0)

    # Test 1: Check if a folium.Map object is returned
    assert isinstance(p, folium.Map)


# plot_landscape
def test_plot_landscape_demo(demo, raster) -> None:
    """Test that the plot_landscape function raises a ValueError if the demography object isn't populated with deme sizes yet"""

    d1 = demography.spDemography()
    with pytest.raises(ValueError):
        plot.plot_landscape(d1, raster, 0)


def test_plot_landscape_timestep(demo, raster) -> None:
    """Test that the plot_landscape function raises a ValueError if the timestep is invalid"""

    with pytest.raises(ValueError):
        plot.plot_landscape(demo, raster, 3)


def test_plot_landscape_returns_plot(demo, raster) -> None:
    """Test that the plot_landscape function returns a matplotlib axes object"""

    p = plot.plot_landscape(demo, raster, 0)
    assert isinstance(p, matplotlib.axes.Axes)


# plot_timeseries
def test_plot_timeseries_data(demo):
    """Test that the plot_timeseries function returns the correct data"""

    times = [0, 1]

    fig, ax = plot.plot_timeseries(demo, times, "Generations")
    x_plot, y_plot = ax.lines[0].get_xydata().T

    # get the y values
    total_individuals = [np.nansum(d) for d in demo.demes]

    assert np.all(x_plot == times)
    assert np.all(y_plot == total_individuals)


def test_plot_timeseries_times_length_1(demo):
    """Test that the plot_timeseries function raises a ValueError if the times list is less than 2"""

    times = [0]

    with pytest.raises(ValueError):
        plot.plot_timeseries(demo, times, "Generations")


def test_plot_timeseries_times_length_not_equal_to_timesteps(demo):
    """Test that the plot_timeseries function raises a ValueError if the length of the times list is not equal to the number of timesteps"""

    times = [0, 1, 2]

    with pytest.raises(ValueError):
        plot.plot_timeseries(demo, times, "Generations")
