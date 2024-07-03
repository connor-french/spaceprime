"""
Module for creating and modifying msprime Demography objects in the service of creating two-dimensional stepping-stone models.
"""

import numpy as np
from typing import Union, List, Optional
import msprime
from .utilities import calc_migration_matrix
import warnings


def add_landscape_change(
    model: msprime.Demography,
    d: np.ndarray,
    timesteps: Union[int, List[int]],
    rate: Union[float, np.ndarray] = 0.001,
    scale: bool = True,
) -> msprime.Demography:
    """
    If there are multiple time steps of deme size change, such as transformed species distribution model projections to past time periods, this function adds the changes to the model.
    It updates the deme sizes and migration rates at each time step.

    Parameters
    ----------
    model : msprime.Demography
        The model object to which the landscape change will be added.
    d : np.ndarray
        The 3D array representing different time steps of population size change.
    timesteps : Union[int, List[int]]
        The list of timesteps representing the amount of time passing between each demographic event, in generations. If a single integer is provided, the function assumes that the time steps are equal.
    rate : Union[float, np.ndarray], optional
        The migration rate. Defaults to 0.001.
    scale : bool, optional
        Whether to scale the migration rate based on population size. Defaults to True.

    Returns
    -------
    msprime.Demography
        The updated model object.
    """

    # check if the number of timesteps matches the number of layers in d, if timesteps is a list
    if isinstance(timesteps, list):
        assert (
            len(timesteps) == d.shape[0] - 1
        ), "Number of timesteps should match the number of layers in d - 1"

    # create a list to store each calculated migration matrix if the user provides a scale for migration rate
    mig_array = []

    # have to set the migration rate differently for the most recent time step, since it's specified prior to add_landscape_change
    if np.array(rate).ndim == 0:
        mig_array.append(calc_migration_matrix(d[0], rate, scale=scale))
    # if the user provides a 3D array of migration matrices, get the migration matrix for the current time step
    else:
        mig_array.append(rate[0])

    # iterate through the first dimension of a 3D array, where the array represents different time steps of population size change
    # omit the most ancient time step (have to set its migration rate differently)
    for step in range(1, d.shape[0]):
        # get the population size values of the current array
        kmat = d[step]

        # get the population size values of array from the more ancient time step
        if step < d.shape[0] - 1:
            kmat_anc = d[step + 1]
        else:
            kmat_anc = None

        # get the population size values array from the more recent time step
        kmat_prev = d[step - 1]

        # get the shape of the array
        n, m = kmat.shape

        # get the time for the demographic event
        if isinstance(timesteps, int):
            demo_time = timesteps * step
        else:
            demo_time = sum(timesteps[:step])

        # if the user provides a scalar migration rate, calculate the migration matrices for all time steps
        if np.array(rate).ndim == 0:
            migmat = calc_migration_matrix(kmat, rate, scale=scale)
        # if the user provides a 3D array of migration matrices, get the migration matrix for the current time step
        else:
            migmat = rate[step]

        # append the migration matrix to the list
        mig_array.append(migmat)

        ##### Update population sizes #####
        # add population size changes according to the values of the current array
        for j in range(n):
            for k in range(m):
                # only update the population size if it is different from the previous time point
                if kmat[j, k] != kmat_prev[j, k]:
                    # add a demographic change to each cell in the raster
                    model.add_population_parameters_change(
                        time=demo_time,
                        population=f"deme_{j}_{k}",
                        initial_size=kmat[j, k],
                    )

        ##### Update migration rates #####
        # add migration rate change for each time step
        # this is updating time steps from the present to the past

        for i in range(n):
            for j in range(m):
                for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                    # check for edges
                    if 0 <= i + di < kmat.shape[0] and 0 <= j + dj < kmat.shape[1]:
                        ## only update migration if the donor and recipient population sizes are different between time points
                        if (
                            kmat_prev[i + di, j + dj] != kmat[i + di, j + dj]
                            and kmat[i, j] != kmat_prev[i, j]
                        ):
                            # get the migration rate from the migration matrix
                            rate = migmat[i + di, j + dj]
                            model.add_migration_rate_change(
                                time=demo_time,
                                rate=rate,
                                source=f"deme_{i}_{j}",
                                dest=f"deme_{i + di}_{j + dj}",
                            )
                        ## have the deme migrate to neighbors if the more ancient time step has an empty deme

                        elif (
                            ## only do this if kmat_anc is not none
                            kmat_anc is not None
                            and kmat_prev[i + di, j + dj] != kmat[i + di, j + dj]
                            and kmat[i, j] != kmat_prev[i, j]
                            and kmat_anc[i, j] <= 1e-9
                        ):
                            # get the migration rate from the migration matrix
                            rate = migmat[i + di, j + dj]
                            model.add_migration_rate_change(
                                time=demo_time,
                                rate=rate,
                                source=f"deme_{i}_{j}",
                                dest=f"deme_{i + di}_{j + dj}",
                            )

    # convert the list of migration matrices to a numpy array
    mig_array = np.array(mig_array)

    # sort the events to make sure everything runs in order
    model.sort_events()

    return model, mig_array


class spDemography(msprime.Demography):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def stepping_stone_2d(
        self,
        d: np.ndarray,
        rate: Union[float, np.ndarray],
        scale: bool = True,
        timesteps: Union[int, List[int]] = None,
    ) -> msprime.Demography:
        """
        Create a 2D stepping stone model, either for a single time step or for multiple time steps of deme size change.

        Parameters
        ----------
        d : numpy.ndarray
            The demography matrix representing the population sizes.
        rate : float or numpy.ndarray
            The migration rate(s) between populations.
            If a float, it represents a constant migration rate for all populations.
            If a numpy.ndarray, it represents a migration matrix with shape (T, N, N),
            where N is the total number of populations and T is the number of time steps - 1, if T > 1.
        scale : bool, optional
            Whether to scale the migration rate matrix. Default is True.
        timesteps : Union[int, List[int]], optional
            The list of timesteps representing the amount of time passing between each demographic event, in generations.
            If a single integer is provided, the function assumes that the time steps are equal.

        Returns
        -------
        model : msprime.Demography
            The constructed 2d stepping stone model as an msprime.Demography object.

        Notes
        -----
        The demography matrix `d` should have shape (n, m) or (k, n, m), where n is the number of rows and m is the number of columns for a 2D array and k is the number of layers in a 3D array.
        The migration rate matrix `rate` should have shape (N, N), where N is the total number of populations.
        If there are multiple time steps of population size change, the `add_landscape_change` function is called to modify the model accordingly.
        """

        if len(d.shape) == 2:
            n = d.shape[0]
            m = d.shape[1]
            N = n * m
            d1 = d  # so I don't have to repeat a bunch of code for the 3D case
        elif len(d.shape) == 3:
            t = d.shape[0]
            n = d.shape[1]
            m = d.shape[2]
            N = n * m
            d1 = d[0]  # so I don't have to repeat a bunch of code for the 3D case

        model = msprime.Demography.isolated_model(d1.reshape(N))

        # set population names
        for j in range(n):
            for k in range(m):
                index = j * m + k
                model.populations[index].name = f"deme_{j}_{k}"

        # setup migration rate matrices and check if the user provided their own migration matrix(es)
        if np.array(rate).ndim == 0:
            if scale:
                model.migration_matrix = calc_migration_matrix(d1, rate, scale=True)
            else:
                model.migration_matrix = calc_migration_matrix(d1, rate, scale=False)
        else:
            if len(rate.shape) == 2:
                assert rate.shape == (
                    N,
                    N,
                ), f"Expected a migration matrix with the shape {(N, N)} and instead got {rate.shape}"
                model.migration_matrix = rate
            else:
                # if the user supplies a 3D array of migration matrices
                assert rate.shape == (
                    t,
                    N,
                    N,
                ), f"Expected a migration matrix with the shape {(t, N, N)} and instead got {rate.shape}"
                model.migration_matrix = rate[0]

        # if there are multiple time steps of population size change
        if len(d.shape) == 3:
            model, mig_array = add_landscape_change(model, d, timesteps, rate, scale)
        else:
            mig_array = model.migration_matrix

        model.sort_events()

        # set the model attributes to the spDemography object
        self.__dict__.update(model.__dict__)
        # return a 3d migration array for each time step
        self.migration_array = mig_array
        # add demes to the object
        self.demes = d

    def add_ancestral_populations(
        self,
        anc_sizes: List[float],
        merge_time: Union[float, int],
        anc_id: Optional[np.ndarray] = None,
        anc_merge_times: Optional[List[float]] = None,
        anc_merge_sizes: Optional[List[float]] = None,
        migration_rate: Optional[float] = None,
    ) -> msprime.Demography:
        """
        Adds ancestral populations to the given demographic model, mapping demes in the spatial simulation to ancestral populations.

        Parameters
        ----------
        anc_sizes : List[float]
            A list of ancestral population sizes.
        merge_time : Union[float, int]
            The time at which all demes in the spatial simulation merge into one or more ancestral populations.
        anc_id : Optional[np.ndarray], optional
            An array of ancestral population IDs- the output of [split_landscape_by_pop][utilities.split_landscape_by_pop]. Defaults to None.
        anc_merge_times : Optional[List[float]], optional
            A list of merge times for ancestral populations. Defaults to None.
        anc_merge_sizes : Optional[List[float]], optional
            A list of sizes for merged ancestral populations. Defaults to None.
        migration_rate : Optional[float], optional
            The symmetric migration rate between ancestral populations. Defaults to None.

        Returns
        -------
        msprime.Demography
            The demographic model with the added ancestral populations.

        Raises
        ------
        ValueError
            If the model already contains ancestral populations.
        ValueError
            If the number of demes in the demographic model does not match the number of demes in the admixture ID raster.

        Notes
        -----
        The function adds ancestral populations to the given demographic model. If `anc_id` is not provided, a single ancestral population is added with the initial size specified in `anc_sizes[0]`. If `anc_id` is provided, a new ancestral population is added for each admixture population, with sizes specified in `anc_sizes`. The demes in the simulation are then merged into their respective ancestral populations based on the values in `anc_id`. If `anc_merge_times` is provided, the ancestral populations are merged at the specified times. If `migration_rate` is provided, symmetric migration is allowed between ancestral populations.
        """

        if any("ANC" in pop.name for pop in self.populations):
            raise ValueError(
                "Model already contains ancestral populations. You need to re-initialize the 2D stepping stone model before adding ancestral populations."
            )

        # convert anc_sizes to a list if it is a single value
        if anc_sizes is not None and not isinstance(anc_sizes, list):
            anc_sizes = [anc_sizes]

        # set anc_id to None and throw a warning if anc_sizes is less than 2
        if len(anc_sizes) < 2 and anc_id is not None:
            warnings.warn(
                "anc_sizes is less than two, but you provided a list of anc_id's. If you meant to setup more than one anc_size, provide a list of sizes to anc_size. Otherwise, you're fine."
            )
            anc_id = None

        # check merge time
        if isinstance(merge_time, list):
            raise ValueError("merge_time should be a single float or int, not a list.")
        elif not isinstance(merge_time, (float, int)):
            raise TypeError("merge_time should be a float or int.")

        # convert anc_merge_times to a list if it is a single value
        if anc_merge_times is not None and not isinstance(anc_merge_times, list):
            anc_merge_times = [anc_merge_times]

        # convert anc_merge_sizes to a list if it is a single value
        if anc_merge_sizes is not None and not isinstance(anc_merge_sizes, list):
            anc_merge_sizes = [anc_merge_sizes]

        if anc_id is None:
            # add an ancestral population
            self.add_population(name="ANC_1", initial_size=anc_sizes[0])

            # get names of populations for initiating the merge
            pop_names = [
                [pop.name] for pop in self.populations if "ANC" not in pop.name
            ]

            # add the time when the spatial simulation collapses into the collecting phase
            [
                self.add_population_split(
                    time=merge_time, derived=name, ancestral="ANC_1"
                )
                for name in pop_names
            ]
        else:
            # make sure the ancestral population ID array has the same number of populations as the array used in the demographic modeling
            try:
                anc_id_1d = anc_id.ravel()
                if len(self.populations) != len(anc_id_1d):
                    raise ValueError
            except ValueError:
                print(
                    f"The number of demes in the demographic model is {len(self.populations)}, while the number of demes in the admixture ID raster is {len(anc_id_1d)}. They need to be the same."
                )

            # add a new ancestral population
            for i in range(1, len(anc_sizes) + 1):
                anc_pop_name = f"ANC_{i}"
                self.add_population(name=anc_pop_name, initial_size=anc_sizes[i - 1])

            # merge each deme in the simulation into its respective ancestral population
            for i in range(len(self.populations)):
                pop_name = self.populations[i].name
                if "ANC" not in pop_name:
                    anc_pop = f"ANC_{int(anc_id_1d[i])}"
                    self.add_population_split(
                        time=merge_time, derived=[pop_name], ancestral=anc_pop
                    )

            # merge ancestral populations at their respective times, if we want that behavior
            if anc_merge_times is not None:
                try:
                    if (
                        len(anc_merge_times) != len(anc_sizes) + 1
                        and len(anc_merge_times) != 1
                    ):
                        raise ValueError
                except ValueError:
                    print("The ancestral merge list should be of length N - 1 or 1")

                # make sure the ancestral merge list is either of size 1 or N - 1
                if len(anc_sizes) > 1 and len(anc_merge_times) > 1:
                    for i in range(1, len(anc_sizes)):
                        if i == 1:
                            n1 = f"ANC_{i}"
                            n2 = f"ANC_{i + 1}"
                            anc_n = f"ANC_{i}_{i + 1}"
                        else:
                            anc_num_str = "_".join(map(str, range(1, i + 1)))
                            n1 = f"ANC_{anc_num_str}"
                            n2 = f"ANC_{i + 1}"
                            anc_n = f"ANC_{anc_num_str}_{i + 1}"
                        ## add the new ancestral populations to the simulation
                        self.add_population(
                            name=anc_n, initial_size=anc_merge_sizes[i - 1]
                        )
                        self.add_population_split(
                            time=anc_merge_times[i - 1],
                            derived=[n1, n2],
                            ancestral=anc_n,
                        )
                ## if there are multiple ancestral populations, but a single merge time
                elif len(anc_sizes) > 1 and len(anc_merge_times) == 1:
                    # get a list of the ancestral population names
                    anc_der_pops = []
                    for i in range(1, len(anc_sizes) + 1):
                        anc_der_name = f"ANC_{i}"
                        anc_der_pops.append(anc_der_name)
                    # make the name of the most ancestral population
                    anc_num_str = "_".join(map(str, range(1, len(anc_sizes) + 1)))
                    anc_n = f"ANC_{anc_num_str}"
                    ## add the new ancestral populations to the simulation
                    self.add_population(name=anc_n, initial_size=anc_merge_sizes[0])
                    self.add_population_split(
                        time=anc_merge_times[0], derived=anc_der_pops, ancestral=anc_n
                    )

            # add migration between ancestral populations
            if migration_rate is not None:
                for i in range(1, len(anc_sizes) + 1):
                    for j in range(i + 1, len(anc_sizes) + 1):
                        self.set_symmetric_migration_rate(
                            populations=[f"ANC_{i}", f"ANC_{j}"],
                            rate=migration_rate,
                        )

        # sort the events to make sure everything runs in order
        self.sort_events()
