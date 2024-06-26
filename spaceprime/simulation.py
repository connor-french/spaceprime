import msprime


def sim_ancestry(*args, **kwargs):
    """
    A thin wrapper around msprime.sim_ancestry. This function takes the same arguments as msprime.simulate
    and calls it directly, allowing users to use simulation functionality within the spaceprime namespace.

    Parameters:
        *args: Variable length argument list.
        **kwargs: Arbitrary keyword arguments.

    Returns:
        The result of msprime.sim_ancestry with the provided arguments.
    """
    return msprime.sim_ancestry(*args, **kwargs)


def sim_mutations(*args, **kwargs):
    """
    A thin wrapper around msprime.sim_mutations. This function takes the same arguments as msprime.simulate
    and calls it directly, allowing users to use simulation functionality within the spaceprime namespace.

    Parameters:
        *args: Variable length argument list.
        **kwargs: Arbitrary keyword arguments.

    Returns:
        The result of msprime.sim_mutations with the provided arguments.
    """
    return msprime.sim_mutations(*args, **kwargs)
