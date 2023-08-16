import numpy as np
import msprime

#### define functions for creating simulations. ####
## NOTE: The add_landscape_change() function is the function that implements demographic changes and needs sped up ##

# create migration matrix from population size matrix or single value
def migration_matrix(populations, rate, scale=True):
    d = populations.shape[0] * populations.shape[1]
    M = np.zeros((d, d))
    
    for i in range(populations.shape[0]):
        for j in range(populations.shape[1]):
            current_index = i * populations.shape[1] + j
            # check the neighboring populations and calculate the migration rates. Looping through tuples was a neat trick! Saved a lot of time.
            for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                # check for edges
                if 0 <= i + di < populations.shape[0] and 0 <= j + dj < populations.shape[1]:
                    neighbor_index = (i + di) * populations.shape[1] + (j + dj)
                    if scale:
                        # mig = donor / recipient * rate unless the pop size is zero
                        M[current_index, neighbor_index] = populations[i + di, j + dj] / populations[i, j] * rate if populations[i, j] != 0 else 0
                    else:
                        M[current_index, neighbor_index] = rate


    return M

# initialize the 2D stepping stone model
def stepping_stone2d(initial_size, rate, scale=True):
    assert len(initial_size.shape) <= 3

    n, m = initial_size.shape
    N = n * m
    model = msprime.Demography.isolated_model(initial_size.reshape(N))

    # set population names
    for j in range(n):
        for k in range(m):
            index = j * m + k
            model.populations[index].name = f"pop_{j}_{k}"

    if np.array(rate).ndim == 0:
        if scale:
            model.migration_matrix = migration_matrix(initial_size, rate, scale=True)
        else: 
            model.migration_matrix = migration_matrix(initial_size, rate, scale=False)
    else:
        assert rate.shape == (N, N), f"Expected a migration matrix with the shape {(N, N)} and instead got {rate.shape}"
        model.migration_matrix = rate

    
    return model

# add demographic changes through time
#### THIS IS THE IMPORTANT FUNCTION ####
def add_landscape_change(model, k_stack, timestep = 1, rate = 0.001, scale=True):
    # iterate through the first dimension of a 3D array, where the array represents different time steps of population size change
    for step in range(1, k_stack.shape[0]):
        # get the population size values of the current array
        kmat = k_stack[step]

        # get the population size values of array from the previous time step
        kmat_prev = k_stack[step - 1]

        # get the shape of the array
        n, m = kmat.shape

        # add population size changes according to the values of the current array
        for j in range(n):
            for k in range(m):
                # only update the population size if it is different from the previous time point
                if kmat[j, k] != kmat_prev[j, k]:
                    # add a demographic change to each cell in the raster
                    model.add_population_parameters_change(time=step * timestep, population=f"pop_{j}_{k}", initial_size=kmat[j, k])

        # add migration rate change for each time step
        ## iterate through the population sizes
        for i in range(n):
            for j in range(m):
                ## also index the neighboring cells
                for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                    ## check for edges
                    if 0 <= i + di < kmat.shape[0] and 0 <= j + dj < kmat.shape[1]:
                        ## only update migration if the donor and recipient population sizes are different between time points
                        if kmat[i + di, j + dj] != kmat_prev[i + di, j +dj] and kmat[i, j] != kmat_prev[i, j]:
                            if scale:
                                ## mig = donor / recipient * rate unless the pop size is zero
                                r = kmat[i + di, j + dj] / kmat[i, j] * rate if kmat[i, j] != 0 else 0
                                model.add_migration_rate_change(time = step * timestep, rate = r, source=f"pop_{i}_{j}", dest=f"pop_{i + di}_{j + dj}")
                            else:
                                model.add_migration_rate_change(time = step * timestep, rate = rate, source=f"pop_{i}_{j}", dest=f"pop_{i + di}_{j + dj}")
                    
                    
    return model



# generate a population size matrix
np.random.seed(2384)
k = np.random.randint(0, 500, (20, 40, 40))


m = stepping_stone2d(k[0], rate = 0.001)

m = add_landscape_change(model = m, k_stack = k, timestep = 1000)

print("Finished!")

