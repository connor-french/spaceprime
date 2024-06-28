"""
Module responsible for pre-processing data for creating demographic models
"""

from os import path
import numpy as np
import numpy.ma as ma
import geopandas as gpd
import pandas as pd
import rasterio
from rasterio.mask import mask
from scipy.interpolate import NearestNDInterpolator
from shapely.geometry import Point
from typing import Union, Optional, List, Tuple, Dict

import allel


## create_raster ##
def create_raster(
    data: np.ndarray,
    reference_raster: rasterio.DatasetReader,
    out_folder: str,
    out_prefix: str,
) -> None:
    """
    Creates a raster dataset from a numpy array and reference raster and writes it to a new GeoTiff file.
    The new raster dataset will have the same dimensions, crs, and transform as the reference raster.

    Parameters
    ----------
    data : np.ndarray
        The numpy array containing the data you want for the raster.
    reference_raster : rasterio.DatasetReader
        The reference rasterio DatasetReader object.
    out_folder : str
        The output folder location where the new raster dataset will be saved.
    out_prefix : str
        The prefix for the output file name.

    Returns
    -------
    None
        The function writes the new raster dataset to the output file location
    """
    with rasterio.open(
        path.join(out_folder, f"{out_prefix}.tif"),
        "w",
        driver="GTiff",
        height=data.shape[0],
        width=data.shape[1],
        count=(
            data.shape[2] if len(data.shape) == 3 else 1
        ),  # check if the data is 2D or 3D
        dtype=data.dtype,
        nodata=-9,
        crs=reference_raster.crs,
        transform=reference_raster.transform,
    ) as dataset:
        if len(data.shape) == 3:  # if data is 3D, write each layer to the dataset
            for i in range(data.shape[2]):
                dataset.write(data[:, :, i], i + 1)
        else:
            dataset.write(data, 1)


def raster_to_demes(
    raster: Union[np.ndarray, rasterio.DatasetReader],
    transformation: str = "linear",
    max_local_size: int = 1000,
    normalize: bool = False,
    threshold: Optional[float] = None,
    thresh_norm: bool = False,
    inflection_point: float = 0.5,
    slope: float = 0.05,
) -> np.ndarray:
    """
    Converts a raster to a 2D np.ndarray of deme sizes using either linear, threshold, or sigmoid transformation functions. For more detail about transformation functions, see [this brief overview](transformation-functions.qmd).
    Raster data should be continuous and positive.
    This function was created with the idea of taking in habitat suitability rasters scaled from 0 to 1, where 0 is no suitability and 1 is the highest suitability.
    However, it is flexible enough to accommodate other continuous rasters that can be coaxed to a 0 to 1 scale with the operation `(data - np.min(data)) / (np.max(data) - np.min(data))` by setting the `normalize` flag to `True`.


    Parameters
    ----------
    raster : Union[np.ndarray, rasterio.DatasetReader]
        The input raster data. It can be a numpy array or a rasterio DatasetReader with one or more layers.
    transformation : str, optional
        The transformation function to be used. Options are "linear", "threshold", and "sigmoid". Default is "linear".
    max_local_size : int, optional
        The maximum local deme size. Default is 1000.
    normalize : bool, optional
        Whether to normalize the raster data. Use if your data is not scaled from 0-1. Default is False.
    threshold : float, optional
        The threshold value for the "threshold" transformation method. Default is None.
    thresh_norm : bool, optional
        Whether to normalize the local deme size based on the average suitability above the threshold. This is useful when comparing thresholded simulations with linear or sigmoid simulations, to maintain similar landscape-wide population sizes across max_local_size values. Default is False.
    inflection_point : float, optional
        The inflection point for the "sigmoid" transformation method. Default is 0.5.
    slope : float, optional
        The slope value for the "sigmoid" transformation method. Default is 0.05.

    Returns
    -------
    np.ndarray
        An ndarray of deme sizes.
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

    # if the input is a rasterio object, read in the raster
    if isinstance(raster, rasterio.DatasetReader):
        d = raster.read(masked=True).filled(0)
    else:
        d = np.where(np.isnan(raster), 0, raster)
    if normalize:

        def _normalize(rast):
            return (rast - np.min(rast)) / (
                np.max(rast) - np.min(rast)
            )  # Normalize the raster data to a 0-1 scale

        d = _normalize(d)  # Apply normalization to the raster data

    if transformation == "linear":
        t = d * max_local_size  # Apply linear transformation to the raster data

        if threshold is not None:
            t[t < max_local_size * threshold] = 1e-10

        t[t > 1e-10] = np.ceil(
            t[t > 1e-10]
        )  # Round up the deme sizes to the nearest integer

    if transformation == "threshold":
        t = d  # Set the deme sizes to the raster data
        if thresh_norm:
            avg_sdm = np.nanmean(
                t[t >= threshold]
            )  # Calculate the average suitability above the threshold
            t[t >= threshold] = np.ceil(
                avg_sdm * max_local_size
            )  # Set the deme sizes above the threshold to the average multiplied by the maximum local size
        else:
            t[t >= threshold] = max_local_size

        t[t < threshold] = 1e-10

    if transformation == "sigmoid":

        def _sigmoid(x, a, b):
            y = 1.0 / (1.0 + np.ma.exp(-(x - a) / b))  # Sigmoid transformation function
            return y

        # vectorize the sigmoid function so it can be applied to the entire raster quickly
        sigmoid_v = np.vectorize(_sigmoid)
        t = (
            sigmoid_v(d, inflection_point, slope) * max_local_size
        )  # Apply sigmoid transformation to the raster data

        if threshold is not None:
            t[t < max_local_size * threshold] = 1e-10

        t[t > 1e-10] = np.ceil(
            t[t > 1e-10]
        )  # Round up the deme sizes to the nearest integer

    t = np.where(t < 1e-10, 1e-10, t)  # Set any deme sizes less than 1e-10 to 1e-10

    return t


######################
## TODO: expand this function take a rasterio DatasetReader object as input for the demes argument####
######################
def calc_migration_matrix(
    demes: np.ndarray, rate: float, scale: bool = True
) -> np.ndarray:
    """
    Calculates a migration matrix based on deme sizes and a global migration rate. The migration rate can be scaled based on population size or set as a constant.

    Parameters
    ----------
    demes : np.ndarray
        The 2D numpy array representing the deme sizes.
    rate : float
        The migration rate.
    scale : bool, optional
        Whether to scale the migration rate based on population size. Default is True.

    Returns
    -------
    np.ndarray
        The migration matrix as a 2D numpy array.
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
                            if (
                                (demes[i + di, j + dj] > 1e-9).all()
                                & (demes[i, j] > 1e-9).all()
                            )
                            else 0
                        )
                    else:
                        # use a constant migration rate
                        M[current_index, neighbor_index] = (
                            rate
                            if (
                                (demes[i + di, j + dj] > 1e-9).all()
                                & (demes[i, j] > 1e-9).all()
                            )
                            else 0
                        )

    return M


## split_landscape_by_pop ##
def split_landscape_by_pop(
    raster: rasterio.DatasetReader,
    coordinates: Union[List[Tuple[float, float]], gpd.GeoDataFrame],
    anc_pop_id: List[Union[int, np.integer]],
    band_index: int = 1,
    mask_rast: bool = False,
) -> np.ma.MaskedArray:
    """
    Uses nearest-neighbor interpolation to classify a landscape raster based on the ancestral population assigned to sampled individuals.
    This function takes in a raster and a list of coordinates and ancestral population IDs assigned to each individual in the empirical data set.
    It then interpolates the population IDs across the landscape and returns the new raster as a masked array.

    Parameters
    ----------
    raster : rasterio.DatasetReader
        The rasterio DatasetReader object representing the landscape raster that you want to divide.
    coordinates : Union[List[Tuple[float, float]], gpd.GeoDataFrame]
        A list of (x, y) coordinates or a geopandas GeoDataFrame representing the coordinates assigned to each individual in the empirical data set.
    anc_pop_id : List[Union[int, np.integer]]
        A list of ancestral population IDs assigned to each empirical individual[^1].
    band_index : int, optional
        The index of the raster to read in. Default is 1. Note- rasterio begins indexing at 1 for raster bands.
    mask_rast : bool, optional
        Whether to mask the interpolation by the landscape. Default is False.

    Returns
    -------
    np.ma.MaskedArray
        The new population assignment raster as a masked array.

    Notes
    -----
    [^1]: These IDs are assigned to each empirical individual typically based on genetic clustering methods like STRUCTURE or PCA. The IDs are used to assign individuals to ancestral populations in the landscape.
    """
    # check if coordinates is a list of tuples or a geopandas GeoDataFrame
    if not isinstance(coordinates, list) and not isinstance(
        coordinates, gpd.GeoDataFrame
    ):
        raise TypeError(
            "The coordinates must be a list of tuples or a geopandas GeoDataFrame."
        )

    if isinstance(coordinates, gpd.GeoDataFrame):
        # convert the geometry column to a list of tuples
        coordinates = list(coordinates.geometry.apply(lambda geom: (geom.x, geom.y)))

    if isinstance(anc_pop_id, pd.Series):
        raise TypeError(
            "The anc_pop_id must be a list of whole numbers. If using a pandas Series, convert it to a list using the to_list() method."
        )

    if not isinstance(anc_pop_id, list) or not all(
        isinstance(x, (int, np.integer)) for x in anc_pop_id
    ):
        raise TypeError("The anc_pop_id must be a list of whole numbers.")

    # read in the raster. This imports the raster as a numpy array
    r2 = raster.read(band_index, masked=True)

    # get the x,y indices of the empirical sampled cells
    indices_x = []
    indices_y = []

    for xy in coordinates:

        # mask requires an iterable as input, so I just repeated the two Point geometries in a list. Mask returns a single value since they overlap in the same place
        pt2 = [Point(xy), Point(xy)]

        # mask the raster with the points
        out_image = mask(raster, pt2, nodata="nan", filled=False)

        # oi first raster
        oi = out_image[0][0]

        # get the locality index
        indices_tup = ma.nonzero(oi)

        indices_x.append(indices_tup[0])
        indices_y.append(indices_tup[1])

    indices_x = np.concatenate(indices_x)
    indices_y = np.concatenate(indices_y)

    # check that the number of indices matches the number of population IDs. If they don't match, this means that some coordinates did not overlap with the raster
    if len(indices_x) != len(anc_pop_id):
        raise ValueError(
            "Some coordinates did not overlap with the raster. Make sure all sampling coordinates overlap with a cell in the raster."
        )

    # get all of the x, y indices of the input raster
    r_x, r_y = np.indices(raster.shape)

    interp = NearestNDInterpolator(list(zip(indices_x, indices_y)), anc_pop_id)
    z = interp(r_x, r_y)
    z = z.astype(int)

    # apply mask of the SDM landscape
    if mask_rast:
        z = ma.masked_array(z, r2.mask, fill_value=-9, dtype=int)

    return z  # Return the new raster as a masked array


def mtp_thresh_from_coords(
    raster: rasterio.DatasetReader,
    coordinates: Union[List[Tuple[float, float]], gpd.GeoDataFrame],
) -> float:
    """
    This function takes the coordinates of empirical sampling localities, finds which raster cells they belong to, extracts the values of the first layer for those localities, and finds the minimum value.
    This value is the maximum threshold value to determine a presence vs absence in a threshold transformation.
    If the threshold is set any higher, empirical sampling localities will not be sampled in the simulations.

    Parameters
    ----------
    raster : rasterio.DatasetReader
        The rasterio DatasetReader object representing the raster data containing the suitability values.
    coordinates : Union[List[Tuple[float, float]], gpd.GeoDataFrame]
        The longitude, latitude coordinates of the empirical sampling localities as a list of coordinate pair tuples or a geopandas GeoDataFrame.

    Returns
    -------
    float
        The maximum threshold value to determine a presence vs absence in a threshold transformation.

    Raises
    ------
    TypeError
        If the coordinates input is not a list, geopandas GeoDataFrame, or pandas DataFrame.
    """
    if isinstance(coordinates, list):
        xy = [Point(xy) for xy in coordinates]
    elif isinstance(coordinates, gpd.GeoDataFrame):
        xy = coordinates.geometry
    else:
        raise TypeError(
            "Invalid coordinates input. Expected list of coordinate pairs or geopandas GeoDataFrame."
        )

    # Mask the first layer with the coordinates
    out_image = mask(raster, xy)
    # Find the minimum value of the masked first layer
    max_thresh = np.nanmin(out_image[0][0])

    return max_thresh


def coords_to_sample_dict(
    raster: Union[np.ndarray, rasterio.DatasetReader],
    coordinates: Union[List[Tuple[float, float]], gpd.GeoDataFrame],
    individual_ids: Optional[List[str]] = None,
    vcf_path: Optional[str] = None,
) -> Tuple[Dict[int, int], Dict[int, np.ndarray], Optional[Dict[int, np.ndarray]]]:
    """
    Convert sample coordinates to sample dictionaries for simulation and analysis. Can optionally include empirical data, which is accepted as a path to a VCF file.

    This function takes a raster, a list of coordinates, and optional individual IDs and VCF path.
    It masks the raster with the given coordinates, retrieves the cell IDs for each individual's locality,
    and returns two dictionaries: a sample dictionary containing the number of individuals to sample from the simulation, and a sample dictionary containing the range of individual indices for each cell ID.
    The first dictionary is used to sample individuals from the simulation, and the second dictionary is used to calculate genetic summary statistics from the sampled individuals.

    Parameters
    ----------
    raster : Union[np.ndarray, rasterio.DatasetReader]
        The raster data as a numpy array or rasterio DatasetReader object.
    coordinates : Union[List[Tuple[float, float]], gpd.GeoDataFrame]
        A list of (x, y) coordinates or a geopandas GeoDataFrame.
    individual_ids : Optional[List[str]], optional
        A list of individual IDs corresponding to those in the VCF file, by default None.
    vcf_path : Optional[str], optional
        The path to the VCF file, by default None.

    Returns
    -------
    Tuple[Dict[int, int], Dict[int, np.ndarray], Optional[Dict[int, np.ndarray]]]
        A tuple containing two or three dictionaries.
        The first dictionary contains the number of individuals to sample from the simulation for each cell ID.
        The second dictionary contains the indices of individuals for each cell ID.
        The third, optional dictionary contains the indices of individuals in the VCF file for each cell ID.

    """
    # check for individual ids and vcf path
    if individual_ids is None:
        try:
            if vcf_path is not None:
                raise ValueError
        except ValueError:
            print(
                "When a VCF path is provided, individual IDs corresponding to those in the VCF are expected. Please provide those IDs."
            )
    elif vcf_path is None:
        try:
            if individual_ids is not None:
                raise ValueError
        except ValueError:
            print(
                "When individual IDs are provided, a VCF file is expected. Please provide a path to your VCF file."
            )
    else:
        deme_dict_empirical = None

    # get the cell that each individual belongs to
    ## I have to iterate through each locality to get the cell ID for each individual locality. Masking everything returns the cell IDs out of order.
    cell_list = []

    # check if coordinates is a geopandas dataframe
    if isinstance(coordinates, gpd.GeoDataFrame):
        # convert the geometry column to a list of tuples
        coordinates = list(coordinates.geometry.apply(lambda geom: (geom.x, geom.y)))

    for xy in coordinates:
        # mask requires an iterable as input, so I just repeated the two Point geometries in a list. Mask returns a single value since they overlap in the same place
        pt2 = [Point(xy), Point(xy)]

        # mask the raster with the points
        out_image = mask(raster, pt2, nodata="nan", filled=False)

        # turn into 1D array
        oi_1d = out_image[0][0].ravel()

        # # get the indices of the locality and append
        cell_id = ma.nonzero(oi_1d)[0]
        cell_list.append(cell_id)

    cell_id_array = np.concatenate(cell_list, axis=0)

    # get the number of individuals to sample from the simulation
    deme_dict_sim = {}
    for cid in np.unique(cell_id_array):
        num_inds = np.sum(cell_id_array == cid)
        deme_dict_sim[int(cid)] = num_inds

    # get the range of indices for each cell id
    deme_dict_sim_long = {}
    unique_ids = np.unique(cell_id_array)
    for i in range(len(unique_ids)):
        cid = unique_ids[i]
        if i > 0:
            cid_prev = unique_ids[i - 1]
            last_ind_prev = deme_dict_sim_long[cid_prev][
                len(deme_dict_sim_long[cid_prev]) - 1
            ]  # get the last index of the previous cell id
            deme_dict_sim_long[int(cid)] = np.array(
                range(last_ind_prev + 1, last_ind_prev + deme_dict_sim[cid] + 1)
            )  # get the range of indices for the current cell id
        else:
            deme_dict_sim_long[int(cid)] = np.array(range(deme_dict_sim[cid]))

    if individual_ids is not None and vcf_path is not None:
        # read in the vcf
        callset = allel.read_vcf(vcf_path)

        # get the indices of each individual in the vcf
        ind_indices = []
        for cid in individual_ids:
            index = list(callset["samples"]).index(cid)
            ind_indices.append(index)

        ind_indices_array = np.array(ind_indices)

        # make each cell id a key and the individuals that belong to that cell id the values
        deme_dict_empirical = {}
        for cid in np.unique(cell_id_array):
            id_inds = np.where(cell_id_array == cid)
            id_ind_indices = ind_indices_array[id_inds]
            deme_dict_empirical[int(cid)] = id_ind_indices
    else:
        deme_dict_empirical = None

    return deme_dict_sim, deme_dict_sim_long, deme_dict_empirical


def anc_to_deme_dict(
    anc_pops: np.ndarray, deme_dict: Dict[int, int]
) -> Dict[int, List[int]]:
    """
    Converts the ancestral population assignments of demes into a dictionary mapping ancestral population IDs to deme indices.

    Parameters
    ----------
    anc_pops : np.ndarray
        An array containing the ancestral population assignments of all demes across the landscape.
    deme_dict : Dict[int, int]
        A dictionary mapping deme indices to the number of individuals being sampled from each deme.

    Returns
    -------
    Dict[int, List[int]]
        A dictionary mapping each ancestral population ID to the range of assigned deme indices.
    """
    # unravel the raster into one dimension for indexing
    r_1d = anc_pops.ravel()

    ap_dict = {}
    for key in deme_dict:
        k = r_1d[key]
        if k in ap_dict.keys():
            ap_dict[k].append(key)
        else:
            ap_dict[k] = [key]

    return ap_dict


def coords_to_deme_dict(
    raster: rasterio.DatasetReader,
    coordinates: Union[List[Tuple[float, float]], gpd.GeoDataFrame],
) -> Dict[int, List[float]]:
    """
    Finds the cells a given set of coordinates belong to in a raster and returns a dictionary mapping the cell indices to the centroid coordinates of those cells.
    Because the cells correspond with demes in the 2D stepping stone models, the cell indices are considered deme indices.
    The coordinates typically correspond to empirical data that the simulations need to be sampled from.

    Parameters
    ----------
    raster : rasterio.DatasetReader
        The raster data as a rasterio DatasetReader object.
    coordinates : Union[List[Tuple[float, float]], gpd.GeoDataFrame]
        A list of (x, y) coordinates or a geopandas GeoDataFrame.

    Returns
    -------
    Dict[int, List[float]]
        A dictionary mapping deme indices to their corresponding coordinates.
    """

    # check if raster is a rasterio.DatasetReader object
    if not isinstance(raster, rasterio.DatasetReader):
        raise TypeError("The raster must be a rasterio.DatasetReader object.")

    # check if coordinates is a list of tuples or a geopandas GeoDataFrame
    if not isinstance(coordinates, list) and not isinstance(
        coordinates, gpd.GeoDataFrame
    ):
        raise TypeError(
            "The coordinates must be a list of tuples or a geopandas GeoDataFrame."
        )

    # check if coordinates is a geopandas dataframe
    if isinstance(coordinates, gpd.GeoDataFrame):
        # convert the geometry column to a list of tuples
        coordinates = list(coordinates.geometry.apply(lambda geom: (geom.x, geom.y)))

    # get the population IDs
    cell_dict = {}
    for xy in coordinates:
        # mask requires an iterable as input, so I just repeated the two Point geometries in a list. Mask returns a single value since they overlap in the same place
        pt2 = [Point(xy), Point(xy)]

        # mask the raster with the points
        out_image = mask(raster, pt2, nodata="nan", filled=False)

        # GET CELL INDICES
        # turn into 1D array
        oi_1d = out_image[0][0].ravel()

        # get the indices of the locality and append
        cell_id = ma.nonzero(oi_1d)[0][0]

        # GET CELL COORDINATES
        # get the indices of the locality
        cell_ind = [tuple(_) for _ in np.transpose(ma.nonzero(out_image[0][0]))]

        rows, cols = zip(
            *cell_ind
        )  # unzip the cell indices to separate row and col lists
        coords_tuple = rasterio.transform.xy(raster.transform, rows, cols)

        # I have to convert the x and y coords to a list instead of a tuple, because the transform.xy function returns a tuple where each long and lat is a single-value list, which makes indexing later convoluted.
        x = coords_tuple[0][0]
        y = coords_tuple[1][0]

        coords = [x, y]

        cell_dict[cell_id] = coords

    return cell_dict
