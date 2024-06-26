import msprime


def sim_ancestry(*args, **kwargs):
    """
    A thin wrapper around msprime.sim_ancestry.

    Args:
      *args: Variable length argument list.
      **kwargs: Arbitrary keyword arguments.

    Returns:
      The TreeSequence object representing the results of the simulation if no replication is performed,
      or an iterator over the independent replicates simulated if the num_replicates parameter has been used.

    Notes:
      This function takes the same arguments as msprime.sim_ancestry and calls it directly,
      allowing users to use simulation functionality within the spaceprime namespace.

      The msprime.sim_ancestry documentation is repeated here for convenience:

      Simulates an ancestral process described by the specified model, demography and
      samples, and return a TreeSequence (or a sequence of replicate tree sequences).

      Parameters:
        samples (Union[int, List[SampleSet], Dict[Union[int, str], int]]): The sampled individuals as either an integer, specifying
          the number of individuals to sample in a single-population model;
          or a list of SampleSet objects defining the properties of
          groups of similar samples; or as a mapping in which the keys
          are population identifiers (either an integer ID or string name)
          and the values are the number of samples to take from the corresponding
          population at its default sampling time. It is important to note that
          samples correspond to individuals here, and each sampled individual
          is usually associated with k sample nodes (or genomes) when
          ploidy = k. See the documentation for more details.
          Either samples or initial_state must be specified.
        demography (Optional[Demography]): The demographic model to simulate, describing the
          extant and ancestral populations, their population sizes and growth
          rates, their migration rates, and demographic events affecting the
          populations over time. See the documentation for more details.
          If not specified (or None) we default to a single population with constant size 1
          (see also the population_size parameter).
        ploidy (int): The number of monoploid genomes per sample individual (Default=2).
          See the documentation for more details.
        sequence_length (Optional[float]): The length of the genome sequence to simulate.
          See the documentation for more details.
        discrete_genome (bool): If True (the default) simulation occurs
          in discrete genome coordinates such that recombination and
          gene conversion breakpoints always occur at integer positions.
          Thus, multiple (e.g.) recombinations can occur at the same
          genome position. If discrete_genome is False simulations
          are performed using continuous genome coordinates. In this
          case multiple events at precisely the same genome location are very
          unlikely (but technically possible).
          See the documentation for more details.
        recombination_rate (Union[float, RateMap]): The rate of recombination along the sequence;
          can be either a single value (specifying a single rate over the entire
          sequence) or an instance of RateMap.
          See the documentation for more details.
        gene_conversion_rate (Optional[float]): The rate of gene conversion along the sequence.
          If provided, a value for gene_conversion_tract_length must also be
          specified. See the documentation for more details.
        gene_conversion_tract_length (Optional[float]): The mean length of the gene conversion
          tracts. For discrete genomes the tract lengths are geometrically
          distributed with mean gene_conversion_tract_length, which must be
          greater than or equal to 1. For continuous genomes the tract lengths are
          exponentially distributed with mean gene_conversion_tract_length,
          which must be larger than 0.
          See the documentation for more details.
        population_size (Optional[int]): The number of individuals of the default single population
          Demography. If not specified, defaults to 1. Cannot be specified
          along with the demography parameter. See the documentation for more details.
        random_seed (Optional[int]): The random seed. If this is not specified or None,
          a high-quality random seed will be automatically generated. Valid random
          seeds must be between 1 and 2^32 - 1.
          See the documentation for more details.
        num_replicates (Optional[int]): The number of replicates of the specified
          parameters to simulate. If this is not specified or None,
          no replication is performed and a TreeSequence object
          returned. If num_replicates is provided, the specified
          number of replicates is performed, and an iterator over the
          resulting TreeSequence objects returned.
          See the documentation for more details.
        record_full_arg (bool): If True, record all intermediate nodes
          arising from common ancestor and recombination events in the output
          tree sequence. This will result in unary nodes (i.e., nodes in marginal
          trees that have only one child). Defaults to False.
          See the documentation for more details.
        additional_nodes (NodeType): Retain all ancestry for any node
          of the specified type. This will result in unary nodes.
          Defaults to class NodeType(0).
          See the documentation for more details.
        coalescing_segments_only (bool): If False, retain all ancestry for any
          nodes that are coalescent nodes anywhere in the sequence. This will result in
          unary nodes. Should be set to False when recording additional nodes.
          Defaults to True.
          See the documentation for more details.
        record_migrations (bool): If True, record all migration events
          that occur in the migration table of
          the output tree sequence. Defaults to False.
          See the documentation for more details.
        initial_state (Optional[TreeSequence]): If specified, initialise the
          simulation from the root segments of this tree sequence and return the
          completed tree sequence. Please see
          the documentation for details of the required
          properties of this tree sequence and its interactions with other parameters.
          All information in the initial_state tables is preserved
          (including metadata) and included in the returned tree sequence.
          (Default: None).
        start_time (Optional[float]): If specified, set the initial time that the
          simulation starts to this value. If not specified, the start
          time is zero if performing a simulation of a set of samples,
          or is the time of the oldest node if simulating from an
          existing tree sequence (see the initial_state parameter).
          See the documentation for more details.
        end_time (Optional[float]): If specified, terminate the simulation at the
          specified time. In the returned tree sequence, all rootward paths from
          samples with time < end_time will end in a node with one child with
          time equal to end_time. Any sample nodes with time >= end_time will
          also be present in the output tree sequence. If not specified or None,
          run the simulation until all samples have an MRCA at all positions in
          the genome. See the documentation for more details.
        record_provenance (bool): If True (the default), record all input
          parameters in the tree sequence provenance.
        model (Union[str, AncestryModel, List[AncestryModel]]): The ancestry model to use. This can be either a
          single instance of AncestryModel (or a string that can be
          interpreted as an ancestry model), or a list of AncestryModel
          instances. If the duration attribute of any of these models is
          set, the simulation will be run until at most t + t_m, where
          t is the simulation time when the model starts and t_m
          is the model's duration. If the duration is not set, the
          simulation will continue until the model completes, the overall
          end_time is reached, or overall coalescence. See
          the documentation for more details,
          and the documentation for the available models
          and examples.

      Returns:
        Union[TreeSequence, Iterator[TreeSequence]]: The TreeSequence object representing the results
        of the simulation if no replication is performed, or an
        iterator over the independent replicates simulated if the
        num_replicates parameter has been used.
    """
    return msprime.sim_ancestry(*args, **kwargs)


def sim_mutations(*args, **kwargs):
    """
    A thin wrapper around msprime.sim_mutations.

    Parameters:
      *args: Variable length argument list.
      **kwargs: Arbitrary keyword arguments.

    Returns:
      The result of msprime.sim_mutations with the provided arguments.

    Notes:
      This function takes the same arguments as msprime.sim_ancestry and calls it directly,
      allowing users to use simulation functionality within the spaceprime namespace.

      The msprime.sim_ancestry parameters are repeated here for convenience:

      Simulates mutations on the specified ancestry and returns the resulting
      tskit.TreeSequence. Mutations are generated at the specified rate
      per unit of sequence length, per unit of time. By default, mutations are
      generated at discrete sites along the genome and multiple mutations
      can occur at any given site. A continuous sequence, infinite-sites model
      can also be specified by setting the `discrete_genome` parameter to
      False.

      If the `model` parameter is specified, this determines the model
      of sequence evolution under which mutations are generated.
      The default mutation model is the msprime.JC69,
      a symmetrical mutation model among the ACGT alleles.
      See the sec_mutations_models section for details of available models.

      If a `random_seed` is specified, this is used to seed the random number
      generator. If the same seed is specified and all other parameters are equal
      then the same mutations will be generated. If no random seed is specified
      then one is generated automatically.

      The time interval over which mutations can occur may be controlled
      using the `start_time` and `end_time` parameters. The `start_time`
      defines the lower bound (in time-ago) on this interval and `max_time`
      the upper bound. Note that we may have mutations associated with
      nodes with time <= `start_time` since mutations store the node at the
      bottom (i.e., towards the leaves) of the branch that they occur on.

      If the tree sequence already has mutations, these are by default retained,
      but can be discarded by passing `keep=False`. However, adding new
      mutations to a tree sequence with existing mutations must be done with
      caution, since it can lead to incorrect or nonsensical results if mutation
      probabilities differ by ancestral state. (As an extreme example, suppose
      that X->Y and X->Z are allowable transitions, but Y->Z is not. If a branch
      already has an X->Y mutation on it, then calling `sim_mutations(...,
      keep=True)` might insert an X->Z mutation above the existing mutation, thus
      implying the impossible chain X->Y->Z.)  However, the effect on nucleotide
      models of mutation are generally very small.

      Note:
        Many mutation models will insert silent transitions (e.g.,
        placing a mutation to A above an existing mutation to A). Such mutations
        are harmless and are required for us to guarantee the statistical
        properties of the process of sequentially adding mutations to a tree
        sequence.

      Parameters:
        tree_sequence (tskit.TreeSequence): The tree sequence we
          wish to throw mutations onto.
        rate (float): The rate of mutation per unit of sequence length
          per unit time, as either a single number (for a uniform rate) or as a
          tskit.RateMap. (Default: 0).
        random_seed (int): The random seed. If this is `None`, a
          random seed will be automatically generated. Valid random
          seeds must be between 1 and 2^32 - 1.
          See the sec_randomness_seeds section for usage examples.
        model (MutationModel): The mutation model to use when generating
          mutations. This can either be a string (e.g., "jc69") or
          an instance of a simulation model class
          e.g, msprime.F84(kappa=0.5).
          If not specified or None, the msprime.JC69
          mutation model is used. Please see the
          sec_mutations_models section for more details
          on specifying mutation models.
        start_time (float): The minimum time ago at which a mutation can
          occur. (Default: no restriction.)
        end_time (float): The maximum time ago at which a mutation can occur
          (Default: no restriction).
        discrete_genome (bool): Whether to generate mutations at only integer positions
          along the genome (Default: True).
        keep (bool): Whether to keep existing mutations. (default: True)

      Returns:
        The tskit.TreeSequence object resulting from overlaying
        mutations on the input tree sequence.
    """
    return msprime.sim_mutations(*args, **kwargs)
