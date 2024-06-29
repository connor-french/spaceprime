"""
Module for post-processing and analyzing simulation results.
"""

import numpy as np
import pandas as pd
import allel
from typing import Tuple, Optional, Dict
import itertools
from math import radians
from sklearn.metrics.pairwise import haversine_distances
from sklearn import linear_model
from esda.moran import Moran
from libpysal.weights import Voronoi
from scipy.stats import entropy


def filter_gt(
    gt: np.ndarray,
    deme_dict_inds: dict = None,
    deme_dict_anc: dict = None,
    missing_data_perc: float = 0,
    r2_thresh: float = 0.1,
    filter_monomorphic: bool = True,
    filter_singletons: bool = True,  # Added new parameter with default value True
) -> Tuple[
    allel.GenotypeArray,
    allel.AlleleCountsArray,
    Optional[Dict[str, allel.AlleleCountsArray]],
    Optional[Dict[str, allel.AlleleCountsArray]],
]:
    """
    Filter genotype matrices output by ts.genotype_matrix() to filter out monomorphic sites, loci in linkage disequilibrium, and recreate missing data patterns common to empirical genotype data.
    Returns the genotype matrix and allele counts matrix for the filtered loci, and optionally allele counts matrices for demes and ancestral populations.

    Parameters
    ----------
    gt : np.ndarray
        The genotype matrix.
    deme_dict_inds : dict, optional
        A dictionary containing the indices of individuals in each deme. Defaults to None.
    deme_dict_anc : dict, optional
        A dictionary containing the indices of individuals in each ancestral population. Defaults to None.
    missing_data_perc : float, optional
        The percentage of missing data allowed. Defaults to 0.
    r2_thresh : float, optional
        The threshold for linkage disequilibrium. Defaults to 0.1.
    filter_monomorphic : bool, optional
        Whether to filter out monomorphic sites, keeping only segregating sites. Defaults to True.
    filter_singletons : bool, optional
        Whether to filter out singletons. Defaults to True.

    Returns
    -------
    Tuple[allel.GenotypeArray, allel.AlleleCountsArray, Optional[Dict[str, allel.AlleleCountsArray]], Optional[Dict[str, allel.AlleleCountsArray]]]
        A tuple containing the filtered genotype matrix, the allele counts matrix, a dictionary of allele counts matrices for demes (if deme_dict_inds is provided), and a dictionary of allele counts matrices for ancestral populations (if deme_dict_anc is provided).

    Notes
    -----
    This function uses a random mask to simulate missing data in the genotype matrix.
    For reproducibility it's advised to set a `np.random.seed()` before calling this function.
    """

    # create data mask that masks the tree sequence genotype matrix to have similar missing data patterns to the empirical data

    # get total the number of elements in the genotype matrix
    total_array_size = gt.shape[0] * gt.shape[1]

    # fill an array with zeros that is the total size of the original matrix
    a = np.zeros(total_array_size, dtype=int)

    # set the first x elements to be the size of non-missing data
    non_miss = int(np.ceil(total_array_size * (1 - missing_data_perc)))

    # 1 means "True"
    a[:non_miss] = 1

    # randomly shuffle True and False values
    np.random.shuffle(a)

    # transform to boolean
    a = a.astype(bool)

    # reshape to genotype matrix shape
    miss_mask = np.reshape(a, gt.shape)

    # mask the genotype matrix
    geno_mat = np.ma.masked_array(gt, mask=miss_mask, fill_value=-1)

    # convert to a GenotypeArray
    h = allel.HaplotypeArray(geno_mat)

    gt = h.to_genotypes(ploidy=2)

    # get the allele counts matrix for all individuals
    ac = gt.count_alleles(max_allele=1)

    if filter_monomorphic:
        # only retain monomorphic alleles
        is_seg = ac.is_segregating()
        gt = gt.compress(is_seg, axis=0)
        ac_seg = ac.compress(is_seg)
    else:
        ac_seg = ac

    # alt allele counts matrix
    # needed for finding unlinked SNPs
    gn = gt.to_n_alt(fill=-1)

    # filter for unlinked loci
    # locate unlinked SNPs
    loc_unlinked = allel.locate_unlinked(gn, threshold=r2_thresh)

    # select unlinked SNPs
    ac_unlinked = ac_seg[loc_unlinked]
    gt_unlinked = gt[loc_unlinked]

    if filter_singletons:
        # remove singletons
        not_singleton = ~ac_unlinked.is_singleton(allele=1)
        ac_unlinked = ac_unlinked[not_singleton]
        gt_unlinked = gt_unlinked[not_singleton]

    if deme_dict_inds is not None:
        ac_unlinked_demes = gt_unlinked.count_alleles_subpops(
            max_allele=1, subpops=deme_dict_inds
        )
    else:
        ac_unlinked_demes = None

    # create allele counts matrices for the ancestral populations
    if deme_dict_anc is not None:
        anc_dict_inds = {}
        for key in deme_dict_anc:
            demes = deme_dict_anc[key]
            anc_dict_inds[key] = []
            for d in demes:
                # convert numpy array to list
                inds = deme_dict_inds[d].tolist()
                anc_dict_inds[key].extend(inds)

        ac_unlinked_anc = gt_unlinked.count_alleles_subpops(
            max_allele=1, subpops=anc_dict_inds
        )

    else:
        ac_unlinked_anc = None

    return gt_unlinked, ac_unlinked, ac_unlinked_demes, ac_unlinked_anc


def calc_sumstats(
    ac: allel.AlleleCountsArray,
    coords_dict: dict,
    anc_demes_dict: dict = None,
    ac_demes: dict = None,
    ac_anc: dict = None,
    between_anc_pop_sumstats: bool = False,
    return_df: bool = False,
    precision: int = 6,
) -> dict:
    """
    Calculates a suite of genetic summary statistics on allele counts matrices generated by `filter_gt` or otherwise generated through scikit-allel.
    The required input is an allele counts matrix for all individuals/demes, and a dictionary mapping deme IDs to their coordinates, generated by `coords_to_deme_dict`.
    Optional inputs are dictionaries mapping ancestral population IDs to their constituent demes (`anc_demes_dict`), an dictionary of allele counts matrices for each deme (`ac_demes`), and dictionary of allele counts matrices for each ancestral population (`ac_anc`).


    Parameters
    ----------
    ac : np.ndarray
        An allele counts matrix for all individuals/demes.
    coords_dict : dict
        A dictionary mapping deme IDs to their coordinates, generated by [coords_to_deme_dict][utilities.coords_to_deme_dict].
    anc_demes_dict : dict, optional
        A dictionary mapping ancestral population IDs to their constituent demes. Defaults to None.
    ac_demes : dict, optional
        A dictionary mapping deme IDs to their allele counts matrices. Necessary if you want to calculate Fst or Dxy between demes. Defaults to None.
    ac_anc : dict, optional
        A dictionary mapping ancestral population IDs to their allel counts matrices. If provided, summary statistics are calculated within ancestral populations and not among them. Defaults to None.
    between_anc_pop_sumstats : bool, optional
        Whether to calculate Fst or Dxy between ancestral populations. Defaults to False.
    return_df : bool, optional
        Whether to return the summary statistics as a pandas DataFrame. Defaults to False.
    precision : int, optional
        The number of decimal places to round the summary statistics to. Defaults to 6.

    Returns
    -------
    dict
        A dictionary of summary statistics.



    Notes
    -----
    This function calculates the following summary statistics, either species-wide or per ancestral population, if provided:
    - Site Frequency Spectrum Hill numbers (q1 and q2), corrected for the number of sites
    - Pi (nucleotide diversity)
    - Tajima's D
    - Pairwise Dxy
        - If `between_anc_pop_sumstats` is True, also calculates pairwise Dxy and Hudson's FST between ancestral populations
    - Pairwise Hudson's FST
        - If `between_anc_pop_sumstats` is True, also calculates pairwise Dxy and Hudson's FST between ancestral populations
    - Isolation-by-distance slope and R2
    - Moran's I

    """

    # calculate scaled Hill numbers of the SFS, using Isaac Overcast's approach
    def _hill_number(freqs, order):
        if order == 0:
            return len(np.nonzero(freqs)[0])
        if order == 1:
            h1 = np.exp(entropy(freqs))
            return h1
        tot = float(np.sum(freqs))
        proportions = np.array(freqs[freqs > 0]) / tot
        prop_order = proportions**order
        h2 = np.sum(prop_order) ** (1 / (1 - order))
        return h2

    # get pairwise dxy for pops
    # deme IDs is a list of population IDs we want to calculate DXY for.
    def _calc_dxy(deme_ids):
        dxy = []
        ## also get the names of each pop pair to identify the raw dxy
        dxy_name = []

        # I'm using the coords_dict dictionary to iterate through the population names, so I can be sure that the DXY distances match up with the geographic distances
        for ac_ind in list(itertools.combinations(deme_ids, 2)):
            ac1 = ac_demes[ac_ind[0]]
            ac2 = ac_demes[ac_ind[1]]

            d = np.nanmean(allel.mean_pairwise_difference_between(ac1, ac2))
            d_name = f"dxy_{ac_ind[0]}_{ac_ind[1]}"

            dxy.append(d)
            dxy_name.append(d_name)
        return dxy, dxy_name

    # isolation-by-distance function
    ##### slope and r2 of gen ~ geo regression ##### (similar to mantel corr)
    def _calc_ibd(dxy, coords):

        # scale dxy according to Roussett 1997 (they used FST, but logic still follows)
        dxy = np.asarray(dxy)

        dxy_scaled = dxy / (1 - dxy)

        # get pairwise geographic distance for all pops

        # convert to radians
        long_rad = [radians(x[0]) for x in coords]
        lat_rad = [radians(y[1]) for y in coords]

        geometry = list(zip(long_rad, lat_rad))

        geometry = [list(_) for _ in geometry]

        dist = (
            haversine_distances(geometry) * 6371000 / 1000
        )  # multiply by Earth radius to get kilometers

        # get the upper diagonal of the distance matrix
        dist = dist[np.triu_indices(dist.shape[0], k=1)]

        # take the log10 of geographic distance to scale linearly for IBD regression
        logdist = np.log10(dist)

        dist_df = pd.DataFrame(data={"geo_dist": logdist, "dxy": dxy_scaled})

        dist_df = dist_df.dropna()

        # linear model, extracting the R2 and coefficient (slope)
        reg = linear_model.LinearRegression()

        # reshape X so it is 2D
        geo_dist = np.array(dist_df["geo_dist"])
        geo_dist = geo_dist.reshape(-1, 1)

        reg.fit(geo_dist, dist_df["dxy"])

        r2 = reg.score(geo_dist, dist_df["dxy"])
        b = reg.coef_[0]

        return b, r2

    # only calculate species-wide stats if there aren't ancestral populations
    if ac_anc is None:
        sfs = allel.sfs_folded(ac)

        # calculate the first 2 Hill numbers of the site frequency spectrum, scaling by the sample size
        sfs_h1 = _hill_number(sfs, 1) / len(sfs)
        sfs_h2 = _hill_number(sfs, 2) / len(sfs)

        # average pi across sites
        pi_raw = allel.mean_pairwise_difference(ac, fill=np.nan)
        pi = np.nanmean(pi_raw)
        # standard deviation of pi across sites
        pi_sd = np.nanstd(pi_raw)

        # tajima's D
        taj_d_raw = allel.moving_tajima_d(ac, size=100, step=10)
        taj_d = np.nanmean(taj_d_raw)
        taj_d_sd = np.nanstd(taj_d_raw)

        stat_dict = {
            "sfs_h1": sfs_h1,
            "sfs_h2": sfs_h2,
            "pi": pi,
            "pi_sd": pi_sd,
            "taj_d": taj_d,
            "taj_d_sd": taj_d_sd,
        }

        if ac_demes is not None:
            # calculate dxy
            dids = list(coords_dict.keys())

            dxy, dxy_name = _calc_dxy(dids)

            for name in range(len(dxy_name)):
                stat_dict[dxy_name[name]] = dxy[name]

    # otherwise, calculate stats per ancestral population
    else:
        stat_dict = {}

        # calc sfs Hill numbers, pi, and tajima's D
        for key in ac_anc:
            sfs = allel.sfs_folded(ac_anc[key])
            stat_dict[f"sfs_h1_a{key}"] = _hill_number(sfs, 1) / len(sfs)
            stat_dict[f"sfs_h2_a{key}"] = _hill_number(sfs, 2) / len(sfs)
            pi_raw = allel.mean_pairwise_difference(ac_anc[key])
            stat_dict[f"pi_a{key}"] = np.nanmean(pi_raw)
            stat_dict[f"pi_sd_a{key}"] = np.nanstd(pi_raw)
            taj_d_raw = allel.moving_tajima_d(ac_anc[key], size=100, step=10)
            stat_dict[f"taj_d_a{key}"] = np.nanmean(taj_d_raw)
            stat_dict[f"taj_d_sd_a{key}"] = np.nanstd(taj_d_raw)

        # Dxy and Hudson's FST for ancestral populations and Dxy for all pops.
        if between_anc_pop_sumstats:
            for ac_ind in list(itertools.combinations(list(ac_anc.keys()), 2)):
                ac1 = ac_anc[ac_ind[0]]
                ac2 = ac_anc[ac_ind[1]]
                d_raw = allel.mean_pairwise_difference_between(ac1, ac2)
                d = np.nanmean(d_raw)
                d_sd = np.nanstd(d_raw)
                stat_dict[f"dxy_a{ac_ind[0]}_{ac_ind[1]}"] = d
                stat_dict[f"dxy_sd_a{ac_ind[0]}_{ac_ind[1]}"] = d_sd
                fst, fst_se, _, _ = allel.average_hudson_fst(ac1, ac2, 100)
                stat_dict[f"fst_a{ac_ind[0]}_{ac_ind[1]}"] = fst
                stat_dict[f"fst_se_a{ac_ind[0]}_{ac_ind[1]}"] = fst_se

            if ac_demes is not None:
                # calculate dxy
                dids = list(coords_dict.keys())

                dxy, dxy_name = _calc_dxy(dids)

                for name in range(len(dxy_name)):
                    stat_dict[dxy_name[name]] = dxy[name]
            else:
                print(
                    "You need to provide allele counts matrices for demes to calculate Dxy between demes."
                )

    # calculate per-deme stats and spatial statistics if ac_demes is provided
    if ac_demes is not None:
        # get pi for all pops
        pi_demes = []
        for key in ac_demes:
            pi_demes.append(np.nanmean(allel.mean_pairwise_difference(ac_demes[key])))

        pi_demes_dict = dict(zip(ac_demes.keys(), pi_demes))

        # add pi to dictionary
        for key in pi_demes_dict:
            colname = f"pi_deme_{key}"
            stat_dict[colname] = pi_demes_dict[key]

            ####### Moran's I ########
        # convert the values of the coords_dict dictionary into a list
        coords = list(coords_dict.values())
        # convert to radians
        long_rad = [radians(x[0]) for x in coords]
        lat_rad = [radians(y[1]) for y in coords]
        ## geographic weights from Voronoi polygons
        coords_rad = np.array([long_rad, lat_rad]).transpose()

        weights = Voronoi(coords_rad, use_index=False)

        # array of per-pop pi, calculated earlier
        pi_array = np.array(pi_demes)
        moran_i = Moran(pi_array, weights)

        mi_i = moran_i.I

        stat_dict["morans_i"] = mi_i

        if anc_demes_dict is None:

            # convert the values of the coords_dict dictionary into a list
            coords = list(coords_dict.values())
            b, r2 = _calc_ibd(dxy, coords)
            stat_dict["ibd_slope"] = b
            stat_dict["ibd_r2"] = r2

        else:
            for key in anc_demes_dict:
                demes = anc_demes_dict[key]

                # get the coordinates associated with the demes
                coords_d = []
                for d in demes:
                    c = coords_dict[d]
                    coords_d.append(c)

                ## calculate dxy. This is super fast, so it's easier to do this than try weird subsetting to subset the full IBD matrix for the DXY values we want per admix pop
                dxy, dxy_name = _calc_dxy(demes)

                for name in range(len(dxy_name)):
                    stat_dict[dxy_name[name]] = dxy[name]

                ## calculate ibd
                b, r2 = _calc_ibd(dxy, coords_d)
                stat_dict[f"ibd_slope_a{key}"] = b
                stat_dict[f"ibd_r2_a{key}"] = r2

    # round all values
    stat_dict = {key: np.round(value, precision) for key, value in stat_dict.items()}

    if return_df:
        stat_df = pd.DataFrame(stat_dict, index=[0])
        return stat_df
    else:

        return stat_dict
