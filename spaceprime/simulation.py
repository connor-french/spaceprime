import msprime


def simulate(*args, **kwargs):
    """
    A thin wrapper around msprime.simulate. This function takes the same arguments as msprime.simulate
    and calls it directly, allowing users to use simulation functionality within the spaceprime namespace.

    Parameters:
        *args: Variable length argument list.
        **kwargs: Arbitrary keyword arguments.

    Returns:
        The result of msprime.simulate with the provided arguments.
    """
    return msprime.simulate(*args, **kwargs)
