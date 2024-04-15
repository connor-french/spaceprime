"""
Module responsible for pre-processing data for creating demographic models
"""

import numpy as np
import rasterio
from typing import Union, Optional
from scipy.interpolate import NearestNDInterpolator
from shapely.geometry import Point
from os.path import splitext
from typing import List, Tuple
from os.path import splitext
import rasterio.mask as mask
import numpy.ma as ma


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

    def check_raster_input(raster):
        if not isinstance(raster, (np.ndarray, rasterio.DatasetReader)):
            raise TypeError(
                "Invalid raster input. Expected numpy array or rasterio DatasetReader object."
            )

    def check_transformation_input(transformation):
        valid_transformations = ["linear", "threshold", "sigmoid"]
        if transformation not in valid_transformations:
            raise ValueError(
                f"Invalid transformation input. Expected one of {valid_transformations}."
            )

    check_raster_input(raster)
    check_transformation_input(transformation)

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


## migration_matrix ##
def migration_matrix(demes: np.ndarray, rate: float, scale: bool = True) -> np.ndarray:
    """
    Calculates a migration matrix based on deme sizes and migration rate. The migration rate can be scaled based on population size or set as a constant.

    **Parameters:**
    - `demes`: The 2D numpy array representing the deme sizes.
    - `rate`: The migration rate.
    - `scale`: Whether to scale the migration rate based on population size. Default is True.

    **Returns:**
    - `M`: The migration matrix as a 2D numpy array.
    """
    if not isinstance(demes, np.ndarray):
        raise TypeError("Invalid demes input. Expected numpy array.")
    if not isinstance(rate, float):
        raise TypeError("Invalid rate input. Expected float.")
    if not isinstance(scale, bool):
        raise TypeError("Invalid scale input. Expected bool.")

    d = demes.shape[0] * demes.shape[1]
    M = np.zeros((d, d))

    for i in range(demes.shape[0]):
        for j in range(demes.shape[1]):
            current_index = i * demes.shape[1] + j
            # check the neighboring demes and calculate the migration rates
            for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                # check for edges
                if 0 <= i + di < demes.shape[0] and 0 <= j + dj < demes.shape[1]:
                    neighbor_index = (i + di) * demes.shape[1] + (j + dj)
                    if scale:
                        # calculate migration rate based on population size
                        # mig = donor / recipient * rate unless the pop size is zero
                        M[current_index, neighbor_index] = (
                            (demes[i + di, j + dj] / demes[i, j]) * rate
                            if demes[i + di, j + dj] > 1e-9 and demes[i, j] > 1e-9
                            else 0
                        )
                    else:
                        # use a constant migration rate
                        M[current_index, neighbor_index] = (
                            rate
                            if demes[i + di, j + dj] > 1e-9 and demes[i, j] > 1e-9
                            else 0
                        )

    return M


## split_landscape_by_pop ##
def split_landscape_by_pop(
    raster_path: str,
    coords: List[Tuple[float, float]],
    admix_id: List[int],
    outfile: str,
    band_index: int = 1,
    mask_rast: bool = False,
) -> None:
    """
    Splits a landscape raster based on the ancestral population assigned to sampled individuals. This function takes in a raster and a list of coordinates and population IDs assigned to each individual in the empirical data set. It then interpolates the population IDs across the landscape and writes the new raster to a file.

    **Parameters:**
    - `raster_path`: The path to the landscape raster that you want to divide.
    - `coords`: A list of tuples representing the coordinates assigned to each individual in the empirical data set.
    - `admix_id`: A list of population IDs assigned to each empirical individual.
    - `outfile`: The path to write the new raster. Must have a .tif extension.
    - `band_index`: The index of the raster to read in. Default is 1. Note- rasterio begins indexing at 1 for raster bands.
    - `mask_rast`: Whether to mask the interpolation by the landscape. Default is False.

    **Returns:**
    - None
    """
    # open the raster
    r = rasterio.open(raster_path)
    # read in the raster. This imports the raster as a numpy array
    r2 = r.read(band_index, masked=True)

    # get the x,y indices of the empirical sampled cells
    indices_x = []
    indices_y = []
    for xy in coords:

        # mask requires an iterable as input, so I just repeated the two Point geometries in a list. Mask returns a single value since they overlap in the same place
        pt2 = [Point(xy), Point(xy)]

        # mask the raster with the points
        out_image = mask(r, pt2, nodata="nan", filled=False)

        # oi first raster
        oi = out_image[0][0]

        # get the locality index
        indices_tup = ma.nonzero(oi)

        indices_x.append(indices_tup[0])
        indices_y.append(indices_tup[1])

    indices_x = np.concatenate(indices_x)
    indices_y = np.concatenate(indices_y)

    # get all of the x, y indices of the input raster
    r_x, r_y = np.indices(r.shape)

    interp = NearestNDInterpolator(list(zip(indices_x, indices_y)), admix_id)
    z = interp(r_x, r_y)

    # apply mask of the SDM landscape
    if mask_rast:
        z = ma.masked_array(z, r2.mask, fill_value=-9, dtype=float)

    # check to make sure the filepath contains the ".tif" suffix, then write the file out
    try:
        root, ext = splitext(outfile)
        if ext != ".tif" or "~" in root:
            raise SyntaxError
        print(f"Landscape with K assignment categories written to {outfile}.")
    except SyntaxError:
        print(
            "The outfile cannot use tilde (~) path expansion and must have a .tif extension."
        )

    with rasterio.open(
        outfile,
        mode="w",
        driver="GTiff",
        height=z.shape[0],
        width=z.shape[1],
        count=1,
        dtype=z.dtype,
        crs=r.crs,
        transform=r.transform,
        nodata=-9,
    ) as new_dataset:
        new_dataset.write(z, 1)
