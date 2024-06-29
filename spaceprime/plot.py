""" Module for visualizing data."""

import numpy as np
import pandas as pd
from shapely.geometry import shape
import geopandas as gpd
import rasterio
import spaceprime
import contextily as ctx
from typing import List
from matplotlib import pyplot as plt


def validate_demographic_model(demo):
    if not hasattr(demo, "demes") or not hasattr(demo, "migration_array"):
        raise ValueError(
            "The demographic model must have 'demes' and 'migration_array' attributes. You may have forgotten to populate the spDemography object with `stepping_stone_2d()`."
        )


def validate_timestep(timestep, demo):
    if timestep < 0 or timestep >= len(demo.demes):
        raise ValueError(f"The timestep must be between 0 and {len(demo.demes) - 1}.")


def get_outgoing_migration_rates(
    demes: np.ndarray, mig_mat: np.ndarray, na_value: float = 1e-9
) -> pd.DataFrame:
    """
    Calculates the outgoing migration rates from each deme based on the migration matrix.

    Parameters
    ----------
    demes : np.ndarray
        An array representing the demes.
    mig_mat : np.ndarray
        A migration matrix representing the migration rates between demes.
    na_value : float, optional
        The value to use for missing migration rates. Defaults to 1e-9.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the outgoing migration rates from each deme.
    """
    # Get the shape of the demes array
    rows, cols = demes.shape
    # Create an empty list to store the migration data
    data = []

    # Iterate over each cell in the demes array
    for i in range(rows):
        for j in range(cols):
            # Calculate the current index based on the row and column
            current_index = i * cols + j
            # Get the deme size at the current cell
            deme_size = demes[i, j]

            # Skip cells with deme size below the specified threshold
            if deme_size <= na_value:
                continue

            # Calculate the migration rates for each direction
            mig_north = mig_mat[current_index, current_index - cols] if i > 0 else 0
            mig_east = mig_mat[current_index, current_index + 1] if j < cols - 1 else 0
            mig_south = (
                mig_mat[current_index, current_index + cols] if i < rows - 1 else 0
            )
            mig_west = mig_mat[current_index, current_index - 1] if j > 0 else 0

            # Append the migration data to the list
            data.append([deme_size, mig_north, mig_east, mig_south, mig_west])

    # Create a pandas DataFrame from the migration data
    return pd.DataFrame(
        data, columns=["deme_size", "mig_north", "mig_east", "mig_south", "mig_west"]
    )


def plot_model(
    demo: spaceprime.spDemography,
    raster: rasterio.DatasetReader,
    timestep: int,
    cmap: str = "viridis",
    legend: bool = True,
    tiles: str = "CartoDB positron",
):
    """
    Plots the demes and migration rates for a given timestep as an interactive map.

    Parameters
    ----------
    demo : spaceprime.spDemography
        The demographic model to plot.
    raster : rasterio.DatasetReader
        The raster dataset used to create the demes matrix(es).
    timestep : int
        The index of the desired timestep to plot.
    cmap : str, optional
        The colormap to use for the deme sizes. Defaults to 'viridis'.
    legend : bool, optional
        Whether to show the color legend. Defaults to True.
    tiles : str, optional
        The basemap tiles to use. Defaults to "CartoDB positron".

    Returns
    -------
    folium.folium.Map
        An interactive plot of the demes and migration rates.

    Example
    -------
    raster = rasterio.open("path/to/raster.tif")
    # Plot the model at timestep 1
    plot_model(demo, raster, 1)

    Notes
    -----
    Since this function returns a `folium` map object, you can further modify the map or save it to an HTML file with the `folium` library.
    """

    # check that the demography object contains demes and migration_array attributes
    validate_demographic_model(demo)

    # check if the number of timesteps is valid
    validate_timestep(timestep, demo)

    # Get the demes matrix for the specified timestep
    demes_matrix = demo.demes[timestep]

    # Mask the array to get only the valid data
    mask = demes_matrix > 1e-10

    # make sure each value is unique so rasterio doesn't merge cells
    array = np.arange(demes_matrix.size).reshape(demes_matrix.shape).astype(np.int32)

    shapes = rasterio.features.shapes(array, mask=mask, transform=raster.transform)

    # read the shapes as separate lists
    dummy_vals = []
    geometry = []
    for shapedict, value in shapes:
        dummy_vals.append(value)
        geometry.append(shape(shapedict))

    # build the gdf object over the two lists
    gdf = gpd.GeoDataFrame({"dummy": dummy_vals, "geometry": geometry}, crs=raster.crs)

    mig_df = get_outgoing_migration_rates(demes_matrix, demo.migration_array[timestep])

    # Add the migration rates to the GeoDataFrame
    gdf = pd.concat([gdf, mig_df], axis=1)

    # Add the deme index to the GeoDataFrame
    gdf["deme_index"] = gdf.index.to_list()

    # Plot the GeoDataFrames
    plot = gdf.explore(
        column="deme_size",
        tooltip=[
            "deme_size",
            "mig_north",
            "mig_east",
            "mig_south",
            "mig_west",
            "deme_index",
        ],
        cmap=cmap,
        legend=legend,
        tiles=tiles,
    )

    return plot


def plot_landscape(
    demo: spaceprime.spDemography,
    raster: rasterio.DatasetReader,
    timestep: int,
    cmap: str = "viridis",
    legend: bool = True,
    basemap: bool = False,
):
    """
    Plots a static map of a transformed landscape at the timestep of your choice.

    Parameters
    ----------
    demo : spaceprime.spDemography
        The demographic model to plot.
    raster : rasterio.DatasetReader
        The raster dataset used to create the demes matrix(es).
    timestep : int
        The timestep to plot.
    cmap : str, optional
        The colormap to use. Defaults to "viridis".
    legend : bool, optional
        Whether to show the colorbar legend. Defaults to True.
    basemap : bool, optional
        Whether to add a basemap. Requires an internet connection. Defaults to False.

    Returns
    -------
    matplotlib.axes.Axes
        A plot of the transformed landscape.

    Note
    ----
    Setting `basemap=True` requires an internet connection to download the basemap tiles. It may take some time to load the tiles depending on your internet speed.
    Since this function returns a `matplotlib` axes object, you can further modify the plot with the `matplotlib` library.
    """

    validate_demographic_model(demo)
    validate_timestep(timestep, demo)

    demes_matrix = demo.demes[timestep]
    mask = demes_matrix > 1e-10

    # plot the demes matrix
    array = np.arange(demes_matrix.size).reshape(demes_matrix.shape).astype(np.int32)

    shapes = rasterio.features.shapes(array, mask=mask, transform=raster.transform)

    dummy_vals = []
    geometry = []
    for shapedict, value in shapes:
        dummy_vals.append(value)
        geometry.append(shape(shapedict))

    gdf = gpd.GeoDataFrame({"dummy": dummy_vals, "geometry": geometry}, crs=raster.crs)
    # add the deme sizes to this gdf, making sure the masked values are omitted
    gdf["deme_size"] = demes_matrix[demes_matrix > 1e-10].flatten()

    if basemap:
        # use gdf.plot to plot the demes matrix
        plot = gdf.plot(column="deme_size", cmap=cmap, legend=legend, alpha=0.7)
        plot = ctx.add_basemap(plot, crs=gdf.crs, source=ctx.providers.CartoDB.Positron)
    else:
        plot = gdf.plot(column="deme_size", cmap=cmap, legend=legend)

    return plot


def plot_timeseries(demo: spaceprime.spDemography, times: List[float], units: str = ""):
    """
    Plots the total number of individuals across the landscape across time.

    Parameters
    ----------
    demo : spaceprime.spDemography
        The demographic model to plot.
    times : List[float]
        A list of times that each landscape timestep corresponds with. This can be in whatever unit you choose.
    units : str, optional
        The units of time the timesteps are specified in. Defaults to a blank string.

    Returns
    -------
    tuple
        A tuple containing the figure and axes objects of the plot.
    """

    validate_demographic_model(demo)

    # check if the times list is greater than 1
    if len(times) < 2:
        raise ValueError("The times list must contain more than one value.")

    # check if the length of the times list is equal to the number of timesteps
    if len(times) != len(demo.demes):
        raise ValueError(
            f"The length of the times list must be equal to the number of timesteps in the demographic model ({len(demo.demes)})."
        )

    # Get the total number of individuals at each timestep
    total_individuals = [np.nansum(d) for d in demo.demes]

    # Plot the timeseries as a line plot, using attractive formatting
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(times, total_individuals, marker="o", color="#224146")
    ax.set_xlabel(f"Time ({units})")
    ax.set_ylabel("Total number of individuals across the landscape")

    return fig, ax
