"""Console script for spaceprime."""

import argparse
import os
import yaml
import rasterio
import msprime
import pandas as pd
import numpy as np
from numpy.random import default_rng
from geopandas import GeoDataFrame
from shapely.geometry import Point
import time
from multiprocessing import Pool

from . import utilities
from . import demography
from . import analysis

# set up logging
import logging

logging.basicConfig(filename="iddc.log", level=logging.INFO)
# write logs to a file


logging.info("spaceprime CLI script started")


# Check if list arguments have more than two elements
def check_argument_length(arg, max_length):
    if arg is not None and len(arg) > max_length:
        raise ValueError(f"{arg} argument can have a maximum of {max_length} elements")


# Check if list arguments are lists, else convert to lists
def check_list_argument(arg):
    if arg is None:
        return None
    elif isinstance(arg, list):
        return arg
    else:
        return [arg]


# Assuming your columns are named 'longitude' and 'latitude'
def read_coords(file_path):
    # check if coords file exists
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Coordinates file {file_path} not found")
    # check if coords file is a CSV
    if not file_path.endswith(".csv"):
        raise ValueError(f"Coordinates file {file_path} must be a CSV")
    df = pd.read_csv(file_path)
    # Check if 'longitude' and 'latitude' columns exist
    if "longitude" not in df.columns or "latitude" not in df.columns:
        raise ValueError("DataFrame must contain 'longitude' and 'latitude' columns.")
    gdf = GeoDataFrame(
        df, geometry=[Point(xy) for xy in zip(df.longitude, df.latitude)]
    )
    return gdf


# read in individual IDs if a path is provided, otherwise use the list
# also perform a series of checks
def read_individuals(args, coords):
    if args.individuals is not None:
        if isinstance(args.individuals, str):
            if not os.path.exists(args.individuals):
                raise FileNotFoundError(
                    f"Individuals file {args.individuals} not found"
                )
            if not args.individuals.endswith(".csv"):
                raise ValueError(f"Individuals file {args.individuals} must be a CSV")
            individuals = pd.read_csv(args.individuals)
            if "individual_id" not in individuals.columns:
                raise ValueError(
                    "Individuals file must contain a column named 'individual_id'"
                )
            individuals = individuals["individual_id"].tolist()
        else:
            individuals = args.individuals

        if len(individuals) != len(coords):
            raise ValueError(
                "Number of individuals must be the same as the number of coordinates"
            )

    return individuals


def read_anc_pop_id(anc_pop_id):
    if anc_pop_id is not None and not len(anc_pop_id) > 1:
        if not os.path.exists(anc_pop_id[0]):
            raise FileNotFoundError(f"anc_pop_id file {anc_pop_id[0]} not found")
        else:
            # Read the CSV file using pandas
            df = pd.read_csv(anc_pop_id[0])
            # Check if 'anc_pop_id' column exists in the CSV file
            if "anc_pop_id" not in df.columns:
                raise ValueError("CSV file does not contain 'anc_pop_id' column")
            # Return the 'anc_pop_id' column as a list
            return df["anc_pop_id"].tolist()
    else:
        return anc_pop_id


# for outputting a map of genetic diversity
def get_map_dict(m, min_num_inds=2, sample_num=2):
    nzp = []

    for i in range(len(m.populations)):
        if (
            m.populations[i].initial_size >= min_num_inds
            and m.populations[i].name != "ANC"
        ):
            nzp.append(m.populations[i].name)

    map_dict = {key: sample_num for key in nzp}
    return map_dict


# get coalescent times for each deme if the map is True
def get_coal_times(tseq, raster, num_anc_pops, sample_num=2, ploidy=2):
    # account for merged ancestral population
    if num_anc_pops > 1:
        num_anc_pops += 1

    coal_list = []
    for j in range(tseq.num_populations):
        # Get the samples corresponding to this population
        samples = tseq.samples(population=j)
        if len(samples) == ploidy * sample_num:
            # Simplify the tree sequence to just these samples
            ts_pop = tseq.simplify(samples=samples)
            div_pi = ts_pop.diversity()
            coal_list.append(div_pi)
        else:
            coal_list.append(-1)

    coal_1d = np.array(coal_list)[:-num_anc_pops]

    coal_array = np.reshape(coal_1d, newshape=raster.shape)

    return coal_array


def setup_demography(
    raster,
    coords,
    max_local_size,
    threshold,
    inflection_point,
    slope,
    mig_rate,
    scale,
    anc_pop_id,
    timesteps,
    anc_sizes,
    merge_time,
    anc_merge_time,
    anc_merge_size,
    anc_mig_rate,
):
    # convert raster to demes
    demes = utilities.raster_to_demes(
        raster=raster,
        max_local_size=max_local_size,
        threshold=threshold,
        inflection_point=inflection_point,
        slope=slope,
    )
    # split landscape by population
    if anc_pop_id is not None:
        anc_pop_mat = utilities.split_landscape_by_pop(
            raster=raster, coordinates=coords, anc_pop_id=anc_pop_id
        )
    else:
        anc_pop_mat = None

    # create demography
    d = demography.stepping_stone_2d(
        d=demes, rate=mig_rate, scale=scale, timesteps=timesteps
    )

    # add ancestral populations
    d = demography.add_ancestral_populations(
        model=d,
        anc_sizes=anc_sizes,
        merge_time=merge_time,
        anc_id=anc_pop_mat,
        anc_merge_times=anc_merge_time,
        anc_merge_sizes=anc_merge_size,
        migration_rate=anc_mig_rate,
    )

    return d


def get_random_value(arg):
    rng = default_rng()
    if arg is None:
        return None
    elif isinstance(arg, (int, float)):
        return arg
    elif len(arg) == 1:
        return arg[0]
    elif isinstance(arg[0], list):
        if isinstance(arg[0][0], int):
            return [rng.integers(a, b + 1) for a, b in arg]
        else:
            return [rng.uniform(a, b) for a, b in arg]
    else:
        if isinstance(arg[0], int):
            return rng.integers(arg[0], arg[1] + 1)
        else:
            return rng.uniform(arg[0], arg[1])


def generate_param_combinations(args):
    param_combos = []
    for _ in range(args.num_param_combos):
        combo = {}
        # max_local_size
        combo["max_local_size"] = get_random_value(args.max_local_size)
        # threshold
        combo["threshold"] = get_random_value(args.threshold)
        # inflection_point
        combo["inflection_point"] = get_random_value(args.inflection_point)
        # slope
        combo["slope"] = get_random_value(args.slope)
        # mig_rate
        combo["mig_rate"] = get_random_value(args.mig_rate)
        # anc_sizes
        if args.anc_sizes is not None:
            combo["anc_sizes"] = get_random_value(args.anc_sizes)
        else:
            combo["anc_sizes"] = None
        # merge_time
        combo["merge_time"] = get_random_value(args.merge_time)
        # anc_merge_time
        if args.anc_merge_time is not None:
            combo["anc_merge_time"] = get_random_value(args.anc_merge_time)
        else:
            combo["anc_merge_time"] = None
        # anc_merge_size
        if args.anc_merge_size is not None:
            combo["anc_merge_size"] = get_random_value(args.anc_merge_size)
        else:
            combo["anc_merge_size"] = None
        # anc_mig_rate
        if args.anc_mig_rate is not None:
            combo["anc_mig_rate"] = get_random_value(args.anc_mig_rate)
        else:
            combo["anc_mig_rate"] = None
        # mutation_rate
        combo["mutation_rate"] = get_random_value(args.mutation_rate)
        # recombination_rate
        combo["recombination_rate"] = get_random_value(args.recombination_rate)

        param_combos.append(combo)

    return param_combos


def run_simulation(combo, args):
    rng = default_rng()
    logging.info(f"Running simulation with parameters: {combo}")
    # read in raster
    r = rasterio.open(args.raster)

    # read in the coordinates
    coords = read_coords(args.coords)

    # read in individuals
    individuals = read_individuals(args, coords)

    # read in anc_pop_id
    anc_pop_id = read_anc_pop_id(args.anc_pop_id)

    # sample dictionaries for ancestry sims and optionally for sumstats
    sample_dicts = utilities.coords_to_sample_dict(r, coords)

    demo_id = f"demo_{rng.integers(0, 2**30)}"
    logging.info(f"Setting up demography with ID: {demo_id}")
    d = setup_demography(
        raster=r,
        coords=coords,
        max_local_size=combo["max_local_size"],
        threshold=combo["threshold"],
        inflection_point=combo["inflection_point"],
        slope=combo["slope"],
        mig_rate=combo["mig_rate"],
        scale=args.scale,
        anc_pop_id=anc_pop_id,
        timesteps=args.timesteps,
        anc_sizes=combo["anc_sizes"],
        merge_time=combo["merge_time"],
        anc_merge_time=combo["anc_merge_time"],
        anc_merge_size=combo["anc_merge_size"],
        anc_mig_rate=combo["anc_mig_rate"],
    )

    print("Finished setting up demography")
    logging.info("Beginning tree sequence simulations")
    start_time = time.time()
    print("Beginning tree sequence simulations")

    # if map is True, return a dictionary mapping samples to all nonzero demes
    # replace the sample dictionary with this dictionary
    if args.map:
        if args.threshold is not None:
            min_num_inds = np.floor(args.threshold * combo["max_local_size"]).astype(
                int
            )
        else:
            min_num_inds = 2

        samples = get_map_dict(d, min_num_inds=min_num_inds)
    else:
        samples = sample_dicts[0]

    for ncoal in range(args.num_coalescent_sims):
        rng = default_rng()
        # set a new random seed for each ancestry simulation
        ancestry_seed = rng.integers(0, 2**30)
        mutation_seed = rng.integers(0, 2**30)

        # simulate tree sequence
        ts = msprime.sim_ancestry(
            samples=samples,
            demography=d,
            sequence_length=args.seq_length,
            recombination_rate=combo["recombination_rate"],
            ploidy=args.ploidy,
            random_seed=ancestry_seed,
            record_provenance=False,
        )

        # simulate mutations
        ts = msprime.sim_mutations(
            ts, rate=combo["mutation_rate"], random_seed=mutation_seed
        )

        end_time = time.time()
        elapsed_time = end_time - start_time

        # write metadata to csv
        metadata = {
            "ancestry_seed": ancestry_seed,
            "mutation_seed": mutation_seed,
            "max_local_size": (
                str(combo["max_local_size"])
                if isinstance(combo["max_local_size"], list)
                else combo["max_local_size"]
            ),
            "threshold": (
                str(combo["threshold"])
                if isinstance(combo["threshold"], list)
                else combo["threshold"]
            ),
            "inflection_point": (
                str(combo["inflection_point"])
                if isinstance(combo["inflection_point"], list)
                else combo["inflection_point"]
            ),
            "slope": (
                str(combo["slope"])
                if isinstance(combo["slope"], list)
                else combo["slope"]
            ),
            "mig_rate": (
                str(combo["mig_rate"])
                if isinstance(combo["mig_rate"], list)
                else combo["mig_rate"]
            ),
            "anc_sizes": (
                str(combo["anc_sizes"])
                if isinstance(combo["anc_sizes"], list)
                else combo["anc_sizes"]
            ),
            "merge_time": (
                str(combo["merge_time"])
                if isinstance(combo["merge_time"], list)
                else combo["merge_time"]
            ),
            "anc_merge_time": (
                str(combo["anc_merge_time"])
                if isinstance(combo["anc_merge_time"], list)
                else combo["anc_merge_time"]
            ),
            "anc_merge_size": (
                str(combo["anc_merge_size"])
                if isinstance(combo["anc_merge_size"], list)
                else combo["anc_merge_size"]
            ),
            "anc_mig_rate": (
                str(combo["anc_mig_rate"])
                if isinstance(combo["anc_mig_rate"], list)
                else combo["anc_mig_rate"]
            ),
            "mutation_rate": (
                str(combo["mutation_rate"])
                if isinstance(combo["mutation_rate"], list)
                else combo["mutation_rate"]
            ),
            "recombination_rate": (
                str(combo["recombination_rate"])
                if isinstance(combo["recombination_rate"], list)
                else combo["recombination_rate"]
            ),
            "demo_id": demo_id,
            "coal_rep": ncoal,
            "sim_time": elapsed_time,
        }
        metadata_file = os.path.join(args.out_folder, f"{args.out_prefix}_metadata.csv")
        if os.path.exists(metadata_file):
            pd.DataFrame(metadata, index=[0]).to_csv(
                metadata_file, mode="a", header=False, index=False
            )
        else:
            pd.DataFrame(metadata, index=[0]).to_csv(metadata_file, index=False)

        logging.info("Beginning coalescent array calculations")
        if args.map:
            if anc_pop_id is not None:
                coal_array = get_coal_times(
                    ts,
                    r,
                    len(set(anc_pop_id)),
                    ploidy=args.ploidy,
                    sample_num=2,
                )
            else:
                coal_array = get_coal_times(ts, r, 1, ploidy=args.ploidy, sample_num=2)

            utilities.create_raster(
                coal_array,
                r,
                out_folder=args.out_folder,
                out_prefix=f"{args.out_prefix}_diversity_map_{ancestry_seed}",
            )
        logging.info("Finished coalescent array calculations and wrote to file")

        # only output other files if args.map is False
        if not args.map:
            # if out_type is 0 or 3, write tree sequence to file
            if args.out_type == 0 or args.out_type == 3:
                ts_file = os.path.join(
                    args.out_folder, f"{args.out_prefix}_ancestry_{ancestry_seed}.trees"
                )
                ts.dump(ts_file)

            # if out_type is 1 or 3, write VCF to file
            if args.out_type == 1 or args.out_type == 3:
                vcf_file = os.path.join(
                    args.out_folder, f"{args.out_prefix}_vcf_{ancestry_seed}.vcf"
                )
                ts.write_vcf(vcf_file, ploidy=args.ploidy, individual_names=individuals)

            # if out_type is 2 or 3, calculate summary statistics
            #### NEED TO FINISH THIS ####
            if args.out_type == 2 or args.out_type == 3:
                # convert tree sequence to a genotype matrix
                gt = ts.genotype_matrix()
                # get necessary dictionaries
                coords_dict = utilities.coords_to_deme_dict(r, coords)

                # split landscape by population if within_anc_pop_sumstats is True or if between_anc_pop_sumstats is True
                if args.within_anc_pop_sumstats or args.between_anc_pop_sumstats:
                    anc_pop_mat = utilities.split_landscape_by_pop(
                        r, coords, anc_pop_id
                    )
                    deme_dict_anc = utilities.anc_to_deme_dict(
                        anc_pop_mat, sample_dicts[1]
                    )
                else:
                    deme_dict_anc = None

                # filter genotype matrix
                _, ac_filt, ac_demes, ac_anc = analysis.filter_gt(
                    gt,
                    deme_dict_inds=sample_dicts[1],
                    deme_dict_anc=deme_dict_anc,
                    missing_data_perc=args.missing_data_perc,
                    r2_thresh=args.r2_thresh,
                    filter_monomorphic=args.filter_monomorphic,
                    filter_singletons=args.filter_singletons,
                )

                # calculate summary statistics
                sumstats = analysis.calc_sumstats(
                    ac=ac_filt,
                    coords_dict=coords_dict,
                    anc_demes_dict=deme_dict_anc,
                    ac_demes=ac_demes,
                    ac_anc=ac_anc,
                    between_anc_pop_sumstats=args.between_anc_pop_sumstats,
                    return_df=True,
                    precision=6,
                )

                # add columns containing the ancestry and mutation seeds
                sumstats["ancestry_seed"] = ancestry_seed
                sumstats["mutation_seed"] = mutation_seed

                # write summary statistics to file
                sumstats_file = os.path.join(
                    args.out_folder, f"{args.out_prefix}_sumstats.csv"
                )
                if os.path.exists(sumstats_file):
                    sumstats.to_csv(sumstats_file, mode="a", header=False, index=False)
                else:
                    sumstats.to_csv(sumstats_file, index=False)

        print("Finished simulating tree sequences")


def main():
    """Console script for spaceprime."""
    parser = argparse.ArgumentParser()
    # global arguments
    parser.add_argument(
        "-p",
        "--params",
        type=str,
        default=None,
        help="Path to YAML file containing parameters. If provided, command line arguments will be ignored. Download a template from the spaceprime documentation. Default is None.",
    )
    parser.add_argument(
        "-r",
        "--raster",
        type=str,
        help="Path to raster file. Can be any raster format that rasterio can read.",
    )
    parser.add_argument(
        "-co",
        "--coords",
        type=str,
        help="Path to file containing sampling coordinates. Must be a CSV file containing columns 'longitude' and 'latitude'.",
    )
    parser.add_argument(
        "-i",
        "--individuals",
        type=lambda x: x.split(",") if "," in x else x,
        default=None,
        help="Either a list or a path to file containing individual IDs. If a CSV, must contain a column named 'individual_id'. Default is None.",
    )
    # Demography Setup
    demography_parser = parser.add_argument_group("Demography Setup")
    demography_parser.add_argument(
        "-n",
        "--normalize",
        type=bool,
        default=False,
        help="Normalize the raster to [0, 1] range. Default is False.",
    )
    demography_parser.add_argument(
        "-t",
        "--transformation",
        type=str,
        default="linear",
        help="Transformation to apply to raster. Options are 'linear', 'threshold', 'sigmoid'. Default is 'linear'.",
    )
    demography_parser.add_argument(
        "-mls",
        "--max_local_size",
        nargs="+",
        type=int,
        default=[1000],
        help="Maximum size of local demes. Accepts a single int or a pair of ints. Default is [1000].",
    )
    demography_parser.add_argument(
        "-th",
        "--threshold",
        nargs="+",
        type=float,
        default=None,
        help="Threshold value for a thresholded transformation. Accepts a single float or a pair of floats. Default is None.",
    )
    demography_parser.add_argument(
        "-ip",
        "--inflection_point",
        nargs="+",
        type=float,
        default=[0.5],
        help="Inflection point value for a sigmoid transformation. Accepts a single int or a pair of ints. Default is [0.5].",
    )
    demography_parser.add_argument(
        "-s",
        "--slope",
        nargs="+",
        type=float,
        default=[0.05],
        help="Slope value for a sigmoid transformation. Accepts a single int or a pair of ints. Default is [0.05].",
    )
    demography_parser.add_argument(
        "-m",
        "--mig_rate",
        nargs="+",
        type=float,
        default=None,
        help="Migration rate between demes. Accepts a single int or a pair of ints. Default is None.",
    )
    demography_parser.add_argument(
        "-sc",
        "--scale",
        type=bool,
        default=True,
        help="Scale the migration rate by donor and recipient deme size using the equation m = (N_donor / N_recipient) * m_global. Default is True.",
    )
    demography_parser.add_argument(
        "-a",
        "--anc_pop_id",
        type=lambda x: x.split(",") if "," in x else x,
        default=None,
        help="List of ancestral population IDs or path to CSV file with column 'anc_pop_id'. Both have to be the same length as --coords. Only supply if you want to add ancestral populations. Default is None.",
    )
    demography_parser.add_argument(
        "-ts",
        "--timesteps",
        type=int,
        default=1,
        help="The list of timesteps representing the amount of time passing between each demographic event, in generations. If a single integer is provided, the function assumes that the time steps are equal. Default is 1.",
    )
    demography_parser.add_argument(
        "-as",
        "--anc_sizes",
        nargs="+",
        type=lambda x: [int(i) for i in x.split(",")],
        default=None,
        help="List of sizes for ancestral populations. Accepts a list of single values or a list of pairs of values. Default is None.",
    )
    demography_parser.add_argument(
        "-mt",
        "--merge_time",
        nargs="+",
        type=int,
        default=None,
        help="Time that demes merge into one or more ancestral populations. Measured in generations. Accepts a single int or a pair of ints. Default is None.",
    )
    demography_parser.add_argument(
        "-amt",
        "--anc_merge_time",
        nargs="+",
        type=int,
        default=None,
        help="Merge time for ancestral populations. Measured in generations. Accepts a single int or a pair of ints. Default is None.",
    )
    demography_parser.add_argument(
        "-ams",
        "--anc_merge_size",
        nargs="+",
        type=int,
        default=None,
        help="Merge size for ancestral populations. Accepts a single int or a pair of ints. Default is None.",
    )
    demography_parser.add_argument(
        "-amr",
        "--anc_mig_rate",
        nargs="+",
        type=float,
        default=None,
        help="Migration rate between ancestral populations. Accepts a single int or a pair of ints. Default is None.",
    )
    # Simulation Setup
    simulation_parser = parser.add_argument_group("Simulation Setup")
    simulation_parser.add_argument(
        "-sl",
        "--seq_length",
        type=int,
        default=1e6,
        help="Length of sequence to simulate. Measured in base pairs. Default is 1e6.",
    )
    simulation_parser.add_argument(
        "-mu",
        "--mutation_rate",
        nargs="+",
        type=float,
        default=[1e-8],
        help="Mutation rate per generation per base pair. Accepts a single float or a pair of floats. Default is [1e-8].",
    )
    simulation_parser.add_argument(
        "-rr",
        "--recombination_rate",
        nargs="+",
        type=float,
        default=[0],
        help="Recombination rate per generation per base pair. Accepts a single float or a pair of floats. Default is [0].",
    )
    simulation_parser.add_argument(
        "-pl",
        "--ploidy",
        type=int,
        default=2,
        help="Ploidy of individuals. Default is 2.",
    )
    simulation_parser.add_argument(
        "-npc",
        "--num_param_combos",
        type=int,
        default=1,
        help="Number of parameter combinations to simulate. Default is 1.",
    )
    simulation_parser.add_argument(
        "-ncs",
        "--num_coalescent_sims",
        type=int,
        default=1,
        help="Number of coalescent simulations to run per parameter combination. Default is 1.",
    )
    # Analysis Setup
    analysis_parser = parser.add_argument_group("Analysis Setup")
    analysis_parser.add_argument(
        "-mdp",
        "--missing_data_perc",
        type=float,
        default=0,
        help="Percentage of data masked as missing. Accepts a single float between 0 and 1. Default is 0.",
    )
    analysis_parser.add_argument(
        "-rt",
        "--r2_thresh",
        type=float,
        default=0.1,
        help="Threshold value for filtering sites based on linkage disequilibrium R^2. Default is 0.1.",
    )
    analysis_parser.add_argument(
        "-fm",
        "--filter_monomorphic",
        type=bool,
        default=True,
        help="Filter out monomorphic sites. Default is True.",
    )
    analysis_parser.add_argument(
        "-fs",
        "--filter_singletons",
        type=bool,
        default=True,
        help="Filter out singleton sites. Default is True.",
    )
    analysis_parser.add_argument(
        "-ss",
        "--sumstats",
        nargs="+",
        type=lambda x: x.split(",") if "," in x else x,
        default="all",
        help="List of summary statistics to calculate. Options are 'pi', 'tajima_d', 'sfs_h', 'fst', 'dxy', 'ibd', 'all'. For explanation of summary statistics, see the documentation. Default is all.",
    )
    analysis_parser.add_argument(
        "-wap",
        "--within_anc_pop_sumstats",
        type=bool,
        default=False,
        help="Whether to calculate summary statistics within ancestral populations. Defaults to False.",
    )
    analysis_parser.add_argument(
        "-bap",
        "--between_anc_pop_sumstats",
        type=bool,
        default=False,
        help="Whether to calculate Fst and/or Dxy between ancestral populations. Defaults to False.",
    )
    # Output
    parser.add_argument(
        "-ot",
        "--out_type",
        type=int,
        default=3,
        help="Type of output. 0 for writing tree sequences to file, 1 for writing VCFs to file, 2 for a CSV of genetic summary statistics, or 3 for writing all outputs to files. Default is 3.",
    )

    parser.add_argument(
        "-map",
        "--map",
        type=bool,
        default=False,
        help="Simulate the genetic diversity for each deme across entire landscape and output a GeoTiff. Overrides out_type. Default is False.",
    )

    parser.add_argument(
        "-of",
        "--out_folder",
        type=str,
        default=None,
        help="Path to output folder. Default is the current working directory.",
    )
    parser.add_argument(
        "-op",
        "--out_prefix",
        type=str,
        default="spaceprime",
        help="Prefix for output files. Default is 'spaceprime'.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        type=bool,
        default=False,
        help="Print progress to console. Default is False.",
    )
    parser.add_argument(
        "-c",
        "--cpu",
        type=int,
        default=1,
        help="Number of CPUs to use for parallel processing. Default is 1.",
    )

    args = parser.parse_args()

    # Check if params YAML file exists. If so, update command line arguments with parameters from YAML file
    if args.params is not None:
        if not os.path.exists(args.params):
            raise FileNotFoundError(f"Parameters file {args.params} not found")
        else:
            with open(args.params, "r") as file:
                params = yaml.safe_load(file)
            # Update command line arguments with parameters from YAML file
            for key, value in params.items():
                setattr(args, key, value)

    args.max_local_size = check_list_argument(args.max_local_size)
    args.inflection_point = check_list_argument(args.inflection_point)
    args.slope = check_list_argument(args.slope)
    args.mig_rate = check_list_argument(args.mig_rate)
    args.anc_pop_id = check_list_argument(args.anc_pop_id)
    args.anc_sizes = check_list_argument(args.anc_sizes)
    # make sure anc_sizes is a list of single value or a list of lists with two values in each element
    if args.anc_sizes[0] is not None:
        if isinstance(args.anc_sizes[0], list):
            args.anc_sizes = [list(map(int, size)) for size in args.anc_sizes]
            for size in args.anc_sizes:
                if len(size) != 2:
                    raise ValueError(
                        f"Each element in --anc_sizes must have two values"
                    )
        else:
            args.anc_sizes = [int(size) for size in args.anc_sizes]

    args.anc_merge_time = check_list_argument(args.anc_merge_time)
    args.anc_merge_size = check_list_argument(args.anc_merge_size)
    args.anc_mig_rate = check_list_argument(args.anc_mig_rate)

    # Check if list arguments that are supposed to be priors have more than two elements
    check_argument_length(args.max_local_size, 2)
    check_argument_length(args.inflection_point, 2)
    check_argument_length(args.slope, 2)
    check_argument_length(args.mig_rate, 2)
    check_argument_length(args.anc_merge_time, 2)
    check_argument_length(args.anc_merge_size, 2)
    check_argument_length(args.anc_mig_rate, 2)

    # check if raster file exists
    if not os.path.exists(args.raster):
        raise FileNotFoundError(f"Raster file {args.raster} not found")

    # check if out_folder exists
    if args.out_folder is not None:
        if not os.path.exists(args.out_folder):
            raise FileNotFoundError(f"Output folder {args.out_folder} not found")

    # generate parameter combinations
    param_combos = generate_param_combinations(args)
    logging.info("Generated parameter combinations")

    # run simulations in parallel
    if args.cpu == 1:
        for combo in param_combos:
            run_simulation(combo, args)
    else:
        with Pool(args.cpu) as p:
            p.starmap(run_simulation, [(combo, args) for combo in param_combos])
