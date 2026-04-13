"""Tests for the `simulation` module."""

import numpy as np
import pytest
import tskit

from spaceprime import demography, simulation


@pytest.fixture
def simple_demo():
    """Minimal spDemography ready to simulate (2x2 grid + single ANC pop)."""
    d = np.array([[1000.0, 500.0], [200.0, 300.0]])
    demo = demography.spDemography()
    demo.stepping_stone_2d(d, rate=0.01)
    demo.add_ancestral_populations(anc_sizes=[5000], merge_time=500)
    return demo


@pytest.fixture
def three_d_demo():
    """Two-timestep spDemography ready to simulate."""
    d = np.array(
        [
            [[1000.0, 500.0], [200.0, 300.0]],
            [[800.0, 400.0], [100.0, 200.0]],
        ]
    )
    demo = demography.spDemography()
    demo.stepping_stone_2d(d, rate=0.01, timesteps=100)
    demo.add_ancestral_populations(anc_sizes=[5000], merge_time=1000)
    return demo


# ---------------------------------------------------------------------------
# sim_ancestry
# ---------------------------------------------------------------------------


class TestSimAncestry:
    """Tests for simulation.sim_ancestry."""

    def test_returns_tree_sequence(self, simple_demo):
        """sim_ancestry returns a tskit.TreeSequence."""
        ts = simulation.sim_ancestry(
            samples={"deme_0_0": 2},
            demography=simple_demo,
            sequence_length=1000,
            random_seed=42,
        )
        assert isinstance(ts, tskit.TreeSequence)

    def test_sample_count(self, simple_demo):
        """Number of individuals in the tree sequence matches the sample dict."""
        samples = {"deme_0_0": 2, "deme_0_1": 3}
        ts = simulation.sim_ancestry(
            samples=samples,
            demography=simple_demo,
            sequence_length=1000,
            random_seed=42,
        )
        assert ts.num_individuals == 5

    def test_multiple_populations_sampled(self, simple_demo):
        """Sampling from multiple demes produces the correct total individual count."""
        samples = {"deme_0_0": 2, "deme_1_0": 2, "deme_1_1": 2}
        ts = simulation.sim_ancestry(
            samples=samples,
            demography=simple_demo,
            sequence_length=5000,
            random_seed=7,
        )
        assert ts.num_individuals == 6

    def test_works_with_3d_model(self, three_d_demo):
        """sim_ancestry succeeds with a multi-timestep demographic model."""
        ts = simulation.sim_ancestry(
            samples={"deme_0_0": 2},
            demography=three_d_demo,
            sequence_length=1000,
            random_seed=1,
        )
        assert isinstance(ts, tskit.TreeSequence)


# ---------------------------------------------------------------------------
# sim_mutations
# ---------------------------------------------------------------------------


class TestSimMutations:
    """Tests for simulation.sim_mutations."""

    def test_returns_tree_sequence(self, simple_demo):
        """sim_mutations returns a tskit.TreeSequence."""
        ts = simulation.sim_ancestry(
            samples={"deme_0_0": 2},
            demography=simple_demo,
            sequence_length=1000,
            random_seed=42,
        )
        mts = simulation.sim_mutations(ts, rate=1e-6, random_seed=42)
        assert isinstance(mts, tskit.TreeSequence)

    def test_nonzero_rate_produces_mutations(self, simple_demo):
        """With a nonzero mutation rate, at least some mutations are produced."""
        ts = simulation.sim_ancestry(
            samples={"deme_0_0": 4, "deme_0_1": 4},
            demography=simple_demo,
            sequence_length=50_000,
            random_seed=42,
        )
        mts = simulation.sim_mutations(ts, rate=1e-6, random_seed=42)
        assert mts.num_mutations > 0

    def test_zero_rate_produces_no_mutations(self, simple_demo):
        """With rate=0, the resulting tree sequence has no mutations."""
        ts = simulation.sim_ancestry(
            samples={"deme_0_0": 2},
            demography=simple_demo,
            sequence_length=1000,
            random_seed=1,
        )
        mts = simulation.sim_mutations(ts, rate=0, random_seed=1)
        assert mts.num_mutations == 0
