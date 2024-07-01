"""
Module providing thin wrappers around msprime.sim_ancestry and msprime.sim_mutations.
"""

import msprime


def sim_ancestry(*args, **kwargs):
    """
    A thin wrapper around msprime.sim_ancestry.

    Parameters
    ----------
    *args : tuple
        Variable length argument list.
    **kwargs : dict
        Arbitrary keyword arguments.

    Returns
    -------
    Union[TreeSequence, Iterator[TreeSequence]]
        The TreeSequence object representing the results of the simulation if no replication is performed,
        or an iterator over the independent replicates simulated if the num_replicates parameter has been used.

    Notes
    -----
    This function takes the same arguments as msprime.sim_ancestry and calls it directly,
    allowing users to use simulation functionality within the spaceprime namespace.

    See the msprime.sim_ancestry documentation for more information: See the msprime documentation for more information: https://tskit.dev/msprime/docs/stable/api.html#msprime.sim_ancestry

    Returns
    -------
    Union[TreeSequence, Iterator[TreeSequence]]
        The TreeSequence object representing the results of the simulation if no replication is performed,
        or an iterator over the independent replicates simulated if the num_replicates parameter has been used.
    """
    return msprime.sim_ancestry(*args, **kwargs)


def sim_mutations(*args, **kwargs):
    """
    A thin wrapper around msprime.sim_mutations.

    Parameters
    ----------
    *args : Variable length argument list.
    **kwargs : Arbitrary keyword arguments.

    Returns
    -------
    The result of msprime.sim_mutations with the provided arguments.

    Notes
    -----
    This function takes the same arguments as msprime.sim_ancestry and calls it directly,
    allowing users to use simulation functionality within the spaceprime namespace.

    See the msprime.sim_mutations documentation for more information: https://tskit.dev/msprime/docs/stable/api.html#msprime.sim_mutations

    """
    return msprime.sim_mutations(*args, **kwargs)
