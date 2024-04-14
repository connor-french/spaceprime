"""
Module responsible for pre-processing data for creating demographic models
"""

import numpy as np
import rasterio
from typing import Union, Optional


## raster_to_demes ##
def raster_to_demes(
    raster: Union[np.ndarray, rasterio.DatasetReader],
    transformation: str = "linear",
    max_local_size: int = 100,
    normalize: bool = False,
    threshold: Optional[float] = None,
    inflection_point: float = 0.5,
    slope: float = 0.05,
) -> np.ndarray:
    """
    Converts a raster to a 2D ndarray of deme sizes using different linear, threshold, or sigmoid transformation functions. For more detail about transformation functions, see [this brief overview](trans-fns.md).
    Raster data should be continuous and positive. This function was created with the idea of taking in habitat suitability rasters scaled from 0 to 1, where 0 is no suitability and 1 is the highest suitability. However, it is flexible enough to accommodate other continuous rasters that can be coaxed to a 0 to 1 scale with the operation `(data - np.min(data)) / (np.max(data) - np.min(data))` by setting the `normalize` flag to `True`.


    **Parameters:**
    - `raster`: The input raster data. It can be a numpy array or a rasterio object with one or more layers.
    - `transformation`: The transformation function to be used. Options are "linear", "threshold", and "sigmoid". Default is "linear".
    - `max_local_size`: The maximum local deme size. Default is 100.
    - `normalize`: Whether to normalize the raster data. Use if your data is not scaled from 0-1. Default is False.
    - `threshold`: The threshold value for the "threshold" transformation method. Default is None.
    - `inflection_point`: The inflection point for the "sigmoid" transformation method. Default is 0.5.
    - `slope`: The slope value for the "sigmoid" transformation method. Default is 0.05.

    **Returns:**
    - `t`: An ndarray of deme sizes.
    """

    d = raster.filled(0)  # Fill any missing values in the raster with 0
    if normalize:

        def _normalize(rast):
            return (rast - np.min(rast)) / (
                np.max(rast) - np.min(rast)
            )  # Normalize the raster data to a 0-1 scale

        d = _normalize(d)  # Apply normalization to the raster data

    if transformation == "linear":
        t = d * max_local_size  # Apply linear transformation to the raster data
        t = np.ceil(t)  # Round up the deme sizes to the nearest integer

    if transformation == "threshold":
        t = d  # Set the deme sizes to the raster data
        avg_sdm = np.nanmean(
            t[t >= threshold]
        )  # Calculate the average suitability above the threshold
        t[t >= threshold] = np.ceil(
            avg_sdm * max_local_size
        )  # Set the deme sizes above the threshold to the average multiplied by the maximum local size
        t[t < threshold] = (
            1e-10  # Set the deme sizes below the threshold to a very small value
        )

    if transformation == "sigmoid":

        def _sigmoid(x, a, b):
            y = 1.0 / (1.0 + np.ma.exp(-(x - a) / b))  # Sigmoid transformation function
            return y

        # vectorize the sigmoid function so it can be applied to the entire raster quickly
        sigmoid_v = np.vectorize(_sigmoid)
        t = (
            sigmoid_v(d, inflection_point, slope) * max_local_size
        )  # Apply sigmoid transformation to the raster data
        t = np.ceil(t)  # Round up the deme sizes to the nearest integer

    t[t == 0] = 1e-10  # Set any deme sizes that are 0 to a very small value

    return t  # Return the ndarray of deme sizes
