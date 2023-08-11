import numpy as np
# to protect divide by zero problems (returns zero if dividing by zero)
def divide_protected(x, y):
    return x / y if y else 0

# create a migration matrix scaled by population size
def k_to_scaled_mig_matrix(k, rate):
    # Convert the input matrix to a 1D array
    k_flat = np.ravel(k)

    # Initialize the migration matrix
    m = np.zeros((k_flat.shape[0], k_flat.shape[0]))

    # Iterate over each cell in the input matrix and compute the migration rates
    for i in range(k_flat.shape[0]):
        for j in range(k_flat.shape[0]):
            if i != j:
                wt_mig = rate * divide_protected(k_flat[j], k_flat[i])
                if (i // 4 == j // 4) and abs(i - j) == 1:
                    # Edge neighbor
                    m[i, j] = wt_mig
                elif (i % 4 == j % 4) and abs(i - j) == 4:
                    # Edge neighbor
                    m[i, j] = wt_mig
                elif abs(i - j) == 3 or abs(i - j) == 12:
                    # Corner neighbor
                    m[i, j] = wt_mig
                elif abs(i - j) == 4:
                    # Center neighbor
                    m[i, j] = wt_mig
    return(m)


# test case 
np.random.seed(10)
krand = np.random.randint(20, size = (10,10))

migmat = k_to_scaled_mig_matrix(krand, 0.001)

migmat