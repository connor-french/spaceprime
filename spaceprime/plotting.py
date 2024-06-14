""" Module for visualizing data."""

import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
import spaceprime as sp


def get_outgoing_migration_rates(
    demes: np.ndarray, mig_mat: np.ndarray, na_value: float = 1e-9
) -> pd.DataFrame:
    """
    Calculates the outgoing migration rates from each deme based on the migration matrix.

    Parameters:
      demes (np.ndarray): An array representing the demes.
      mig_mat (np.ndarray): A migration matrix representing the migration rates between demes.
      na_value (float, optional): The value to use for missing migration rates. Defaults to 1e-9.

    Returns:
      pd.DataFrame: A DataFrame containing the outgoing migration rates from each deme.
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
    demo: sp.spDemography,
    raster: rasterio.DatasetReader,
    timestep: int,
    cmap: str = "viridis",
    legend: bool = True,
    tiles: str = "CartoDB positron",
):
    """
    Plots the demes and migration rates for a given timestep as an interactive map.

    Parameters:
      demo (sp.spDemography): The demographic model to plot.
      raster (rasterio.DatasetReader): The raster dataset used to create the demes matrix(es).
      timestep (int): The index of the desired timestep to plot.
      cmap (str, optional): The colormap to use for the deme sizes. Defaults to 'viridis'.
      legend (bool, optional): Whether to show the color legend. Defaults to True.
      tiles (str, optional): The basemap tiles to use. Defaults to "CartoDB positron".

    Returns:
      folium.folium.Map: An interactive plot of the demes and migration rates.

    Example:
      raster = rasterio.open("path/to/raster.tif")
      # Plot the model at timestep 1
      plot_model(demo, raster, 1)

    """

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
