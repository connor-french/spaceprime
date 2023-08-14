import numpy as np

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
                        # in case someone doesn't want a scaled rate
                        M[current_index, neighbor_index] = rate


    return M


# test case 
np.random.seed(10)
krand = np.random.randint(20, size = (5,5))

migmat = migration_matrix(krand, 0.001)

migmat