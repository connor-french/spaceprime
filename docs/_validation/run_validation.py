#!/usr/bin/env python3
"""Validate spaceprime simulations against stepping-stone population-genetic theory.

This script builds demographic models *with spaceprime* (the thing under test) and
analyses the resulting coalescent simulations with independent, transparent code
(scikit-allel + numpy) so the validation does not lean on the same package's
summary-statistic functions.

Experiments
-----------
A. Dimensionality signature of isolation by distance (Rousset 1997):
   linearized FST should be ~linear in geographic distance in a 1D habitat and
   ~linear in log(distance) in a 2D habitat.
B. Neighborhood size controls the IBD slope (Rousset 1997): increasing Nm should
   monotonically decrease both the IBD slope and the mean FST in a 1D chain.
C. Heterogeneous landscape / isolation by resistance (McRae 2006): in a gridded
   landscape with a low-suitability barrier, genetic differentiation should be
   better predicted by resistance distance (from the migration graph) than by
   straight-line Euclidean distance.

Outputs: PNG figures in figures/ and a machine-readable results.json next to this
file. Run under the dev env from the repo root:
    pixi run -e dev python docs/_validation/run_validation.py
"""

from __future__ import annotations

import json
import time
from itertools import combinations
from pathlib import Path

import allel
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import spaceprime as sp
from spaceprime import demography, simulation, utilities

HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
FIG.mkdir(parents=True, exist_ok=True)

RNG_SEED = 20240611
MUT_RATE = 1e-7
RECOMB_RATE = 1e-8


# --------------------------------------------------------------------------- #
# Simulation + analysis helpers
# --------------------------------------------------------------------------- #
def build_model(deme_sizes: np.ndarray, rate: float, anc_size: float, merge_time: float):
    """Build a spaceprime 2D stepping-stone model with a deep ancestral backstop.

    The ancestral merge only caps the (rare) deepest lineages so every run
    terminates quickly; it is far older than the within-landscape structure of
    interest and does not affect the relative differentiation among demes.
    """
    demo = demography.spDemography()
    demo.stepping_stone_2d(deme_sizes, rate=rate, scale=True)
    demo.add_ancestral_populations(anc_sizes=[anc_size], merge_time=merge_time)
    return demo


def occupied_cells(deme_sizes: np.ndarray, threshold: float = 1e-9):
    """Return (row, col) of demes large enough to sample/inhabit."""
    return [
        (i, j)
        for i in range(deme_sizes.shape[0])
        for j in range(deme_sizes.shape[1])
        if deme_sizes[i, j] > threshold
    ]


def simulate(deme_sizes, rate, sample_cells, n_diploid, seq_len, anc_size,
             merge_time, seed):
    """Simulate a tree sequence; return (mts, ordered list of sampled cells)."""
    demo = build_model(deme_sizes, rate, anc_size, merge_time)
    samples = {f"deme_{i}_{j}": n_diploid for (i, j) in sample_cells}
    ts = simulation.sim_ancestry(
        samples=samples,
        demography=demo,
        sequence_length=seq_len,
        recombination_rate=RECOMB_RATE,
        random_seed=seed,
    )
    mts = simulation.sim_mutations(ts, rate=MUT_RATE, random_seed=seed + 1)
    return mts, demo


def allele_counts_per_deme(mts, deme_sizes):
    """Map each sampled deme -> scikit-allel AlleleCountsArray over segregating sites.

    Sample nodes are grouped by their msprime population id, then the population
    id is decoded back to its (row, col) via the population name 'deme_i_j'.
    """
    pop_name = {p.id: p.metadata.get("name", "") if p.metadata else ""
                for p in mts.populations()}
    # msprime stores the name on the population table; fall back to that.
    names = [mts.population(p).metadata.get("name", "") for p in range(mts.num_populations)]

    samples = mts.samples()
    samp_pop = np.array([mts.node(u).population for u in samples])
    gm = mts.genotype_matrix()  # (n_sites, n_samples), haplotype (0/1)

    ac_by_cell = {}
    for pid in np.unique(samp_pop):
        nm = names[pid]
        if not nm.startswith("deme_"):
            continue
        _, i, j = nm.split("_")
        cols = np.where(samp_pop == pid)[0]
        h = allel.HaplotypeArray(gm[:, cols])
        ac_by_cell[(int(i), int(j))] = h.count_alleles(max_allele=1)
    return ac_by_cell


def pairwise_table(ac_by_cell, cells):
    """Pairwise linearized Hudson FST and dxy for every pair of sampled cells."""
    rows = []
    for a, b in combinations(cells, 2):
        ac1, ac2 = ac_by_cell[a], ac_by_cell[b]
        fst, _, _, _ = allel.average_hudson_fst(ac1, ac2, blen=100)
        fst = max(fst, 0.0)
        dxy = np.nanmean(allel.mean_pairwise_difference_between(ac1, ac2))
        rows.append(dict(a=a, b=b, fst=fst, fst_lin=fst / (1 - fst) if fst < 1 else np.nan,
                         dxy=dxy))
    return rows


def lattice_dist(a, b):
    return np.hypot(a[0] - b[0], a[1] - b[1])


def linregress(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    slope, intercept = np.polyfit(x, y, 1)
    yhat = slope * x + intercept
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return slope, intercept, r2


# --------------------------------------------------------------------------- #
# Experiment A: 1D vs 2D dimensionality signature (Rousset 1997)
# --------------------------------------------------------------------------- #
def averaged_pairs(deme_sizes, rate, cells, n_diploid, seq_len, seed0, reps):
    """Run `reps` independent simulations and average per-pair fst_lin and dxy.

    Pairwise Hudson FST from a handful of diploids per deme is noisy; averaging
    the same pair across independent replicates sharpens the distance signal so
    the *shape* (linear vs log) can be resolved, mirroring the multi-replicate
    approach of Szep et al. (2022).
    """
    acc = {}
    for r in range(reps):
        mts, _ = simulate(deme_sizes, rate, cells, n_diploid, seq_len,
                          anc_size=1000, merge_time=200000, seed=seed0 + 17 * r)
        ac = allele_counts_per_deme(mts, deme_sizes)
        for row in pairwise_table(ac, cells):
            key = (row["a"], row["b"])
            acc.setdefault(key, {"fst_lin": [], "dxy": []})
            acc[key]["fst_lin"].append(row["fst_lin"])
            acc[key]["dxy"].append(row["dxy"])
    pairs = sorted(acc.keys())
    dist = np.array([lattice_dist(a, b) for (a, b) in pairs])
    flin = np.array([np.nanmean(acc[k]["fst_lin"]) for k in pairs])
    dxy = np.array([np.nanmean(acc[k]["dxy"]) for k in pairs])
    return dist, flin, dxy


def experiment_a():
    print("\n=== Experiment A: dimensionality signature (Rousset 1997) ===")
    N = 300
    rate = 2.0 / N  # Nm ~ 2 per neighbour
    reps = 3
    out = {}

    # ---- 1D linear habitat: 1 x 30 ----
    d1 = np.full((1, 30), float(N))
    dist1, flin1, dxy1 = averaged_pairs(d1, rate, occupied_cells(d1),
                                        n_diploid=8, seq_len=5e6,
                                        seed0=RNG_SEED, reps=reps)
    # ---- 2D habitat: 10 x 10 ----
    d2 = np.full((10, 10), float(N))
    dist2, flin2, dxy2 = averaged_pairs(d2, rate, occupied_cells(d2),
                                        n_diploid=6, seq_len=5e6,
                                        seed0=RNG_SEED + 1000, reps=reps)

    for tag, dist, flin, dxy, n in [("1d", dist1, flin1, dxy1, 30),
                                    ("2d", dist2, flin2, dxy2, 100)]:
        _, _, r2_f_d = linregress(dist, flin)
        _, _, r2_f_l = linregress(np.log(dist), flin)
        _, _, r2_x_d = linregress(dist, dxy)
        _, _, r2_x_l = linregress(np.log(dist), dxy)
        out[tag] = dict(n_demes=n, n_pairs=int(len(dist)), reps=reps,
                        fst_r2_vs_dist=r2_f_d, fst_r2_vs_logdist=r2_f_l,
                        dxy_r2_vs_dist=r2_x_d, dxy_r2_vs_logdist=r2_x_l)
        print(f"  {tag.upper()}: FST/(1-FST)  linear R2={r2_f_d:.3f}  log R2={r2_f_l:.3f}"
              f"  |  dxy  linear R2={r2_x_d:.3f}  log R2={r2_x_l:.3f}")

    # ---- figure: FST panels on top row, 2D dxy centered on bottom row ----
    fig = plt.figure(figsize=(14, 9))
    gs = fig.add_gridspec(2, 4, hspace=0.45, wspace=0.35)
    ax_1d  = fig.add_subplot(gs[0, 0:2])   # 1D FST (top-left half)
    ax_2d  = fig.add_subplot(gs[0, 2:4])   # 2D FST (top-right half)
    ax_dxy = fig.add_subplot(gs[1, 1:3])   # 2D dxy (bottom centre)

    for ax, dist, y, ylabel, title in [
        (ax_1d,  dist1, flin1,
         r"$F_{ST}/(1-F_{ST})$  (mean of 3 reps)", "1D habitat (1 × 30)"),
        (ax_2d,  dist2, flin2,
         r"$F_{ST}/(1-F_{ST})$  (mean of 3 reps)", "2D habitat (10 × 10)"),
        (ax_dxy, dist2, dxy2,
         r"$d_{xy}$  (mean of 3 reps)",             "2D habitat (10 × 10)"),
    ]:
        ax.scatter(dist, y, s=10, alpha=0.30, color="#3b6", edgecolor="none")
        xs = np.linspace(dist.min(), dist.max(), 100)
        sl_d, ic_d, r2d = linregress(dist, y)
        ax.plot(xs, sl_d * xs + ic_d, "b-", lw=1.7,
                label=f"linear in dist  R$^2$={r2d:.2f}")
        sl_l, ic_l, r2l = linregress(np.log(dist), y)
        ax.plot(xs, sl_l * np.log(xs) + ic_l, "r--", lw=1.7,
                label=f"linear in log(dist)  R$^2$={r2l:.2f}")
        ax.set_xlabel("geographic distance (deme steps)")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(frameon=False, fontsize=9)

    fig.suptitle("Isolation by distance: 1D is linear in distance, 2D in log-distance "
                 "(Rousset 1997)", fontsize=12)
    fig.savefig(FIG / "exp_a_dimensionality.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("  wrote figures/exp_a_dimensionality.png")
    return out


# --------------------------------------------------------------------------- #
# Experiment B: neighborhood size controls IBD slope (Rousset 1997)
# --------------------------------------------------------------------------- #
def experiment_b():
    print("\n=== Experiment B: IBD slope vs neighborhood size (Nm) ===")
    N = 200
    n_dip = 6
    nm_values = [0.5, 1.0, 2.0, 4.0, 8.0]
    d = np.full((1, 25), float(N))
    cells = occupied_cells(d)
    results = []
    for k, nm in enumerate(nm_values):
        rate = nm / N
        mts, _ = simulate(d, rate, cells, n_dip, seq_len=3e6, anc_size=1000,
                          merge_time=200000, seed=RNG_SEED + 200 + k)
        ac = allele_counts_per_deme(mts, d)
        tab = pairwise_table(ac, cells)
        dist = np.array([lattice_dist(r["a"], r["b"]) for r in tab])
        flin = np.array([r["fst_lin"] for r in tab])
        slope, _, r2 = linregress(dist, flin)
        mean_fst = float(np.nanmean([r["fst"] for r in tab]))
        nbr = [r["fst"] for r in tab if abs(lattice_dist(r["a"], r["b"]) - 1) < 1e-9]
        results.append(dict(nm=nm, rate=rate, ibd_slope=slope, ibd_r2=r2,
                            mean_fst=mean_fst, nbr_fst=float(np.nanmean(nbr))))
        print(f"  Nm={nm:>4}: slope={slope:.4e}  meanFST={mean_fst:.3f}  R2={r2:.2f}")

    # ---- figure ----
    nm = [r["nm"] for r in results]
    slopes = [r["ibd_slope"] for r in results]
    fsts = [r["mean_fst"] for r in results]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))
    axes[0].plot(nm, slopes, "o-", color="#b34")
    axes[0].set_xlabel("neighborhood size  Nm (per neighbour)")
    axes[0].set_ylabel("IBD slope  d[$F/(1-F)$]/d(distance)")
    axes[0].set_title("IBD slope decreases with gene flow")
    axes[0].set_xscale("log")
    axes[1].plot(nm, fsts, "s-", color="#36b")
    axes[1].set_xlabel("neighborhood size  Nm (per neighbour)")
    axes[1].set_ylabel(r"mean pairwise $F_{ST}$")
    axes[1].set_title("Differentiation decreases with gene flow")
    axes[1].set_xscale("log")
    fig.tight_layout()
    fig.savefig(FIG / "exp_b_slope_vs_nm.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("  wrote figures/exp_b_slope_vs_nm.png")
    return dict(results=results)


# --------------------------------------------------------------------------- #
# Experiment C: heterogeneous landscape / isolation by resistance (McRae 2006)
# --------------------------------------------------------------------------- #
def resistance_distances(deme_sizes, rate, cells):
    """Effective resistance between sampled cells from the migration graph.

    Conductance between neighbours = symmetric migration rate (spaceprime's
    size-scaled rate). Resistance distance R_ij = L+_ii + L+_jj - 2 L+_ij on
    the connected component's Laplacian pseudoinverse.
    """
    M = utilities.calc_migration_matrix(deme_sizes, rate, scale=True)
    C = 0.5 * (M + M.T)  # symmetric conductance
    deg = C.sum(axis=1)
    nodes = np.where(deg > 0)[0]  # drop isolated (barrier) demes
    idx = {n: k for k, n in enumerate(nodes)}
    Csub = C[np.ix_(nodes, nodes)]
    L = np.diag(Csub.sum(axis=1)) - Csub
    Lp = np.linalg.pinv(L)
    m = deme_sizes.shape[1]
    out = {}
    for a, b in combinations(cells, 2):
        ia = idx[a[0] * m + a[1]]
        ib = idx[b[0] * m + b[1]]
        out[(a, b)] = Lp[ia, ia] + Lp[ib, ib] - 2 * Lp[ia, ib]
    return out


def experiment_c():
    print("\n=== Experiment C: isolation by resistance (McRae 2006) ===")
    # Build a suitability raster: suitable everywhere except a central vertical
    # barrier column, with a one-cell corridor so the landscape stays connected.
    nrow, ncol = 11, 11
    suit = np.ones((nrow, ncol), dtype=float)
    barrier_col = 5
    suit[:, barrier_col] = 0.0
    corridor_row = nrow // 2
    suit[corridor_row, barrier_col] = 1.0  # gap in the barrier

    # Use spaceprime's raster->deme transform (linear) to get deme sizes.
    deme_sizes = utilities.raster_to_demes(suit, transformation="linear",
                                           max_local_size=200)
    cells = occupied_cells(deme_sizes)
    rate = 2.0 / 200  # Nm ~ 2 on suitable cells

    # Sample a manageable subset: every cell two columns from the barrier on both
    # sides plus the edges, across all rows -> good cross-barrier coverage.
    sample_cols = [0, 2, 3, 7, 8, 10]
    sample_cells = [(i, j) for (i, j) in cells if j in sample_cols]

    mts, _ = simulate(deme_sizes, rate, sample_cells, n_diploid=6, seq_len=4e6,
                      anc_size=1000, merge_time=300000, seed=RNG_SEED + 300)
    ac = allele_counts_per_deme(mts, deme_sizes)
    tab = pairwise_table(ac, sample_cells)

    rdist = resistance_distances(deme_sizes, rate, sample_cells)
    eucl = np.array([lattice_dist(r["a"], r["b"]) for r in tab])
    res = np.array([rdist[(r["a"], r["b"])] for r in tab])
    flin = np.array([r["fst_lin"] for r in tab])
    cross = np.array([(r["a"][1] < barrier_col) != (r["b"][1] < barrier_col)
                      for r in tab])

    _, _, r2_eucl = linregress(eucl, flin)
    _, _, r2_res = linregress(res, flin)
    print(f"  R2(genetic ~ Euclidean)   = {r2_eucl:.3f}")
    print(f"  R2(genetic ~ resistance)  = {r2_res:.3f}")
    print(f"  cross-barrier pairs: {cross.sum()} / {len(cross)}")

    # ---- figures: landscape + the two scatterplots ----
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    im = axes[0].imshow(np.where(deme_sizes > 1e-9, deme_sizes, np.nan),
                        cmap="YlGn", origin="upper")
    sc_r = [c[0] for c in sample_cells]
    sc_c = [c[1] for c in sample_cells]
    axes[0].scatter(sc_c, sc_r, c="k", s=18, marker="x")
    axes[0].set_title("Landscape (deme size) + sampled cells")
    fig.colorbar(im, ax=axes[0], shrink=0.8, label="deme size")

    for ax, x, lab, r2 in [
        (axes[1], eucl, "Euclidean distance", r2_eucl),
        (axes[2], res, "resistance distance", r2_res),
    ]:
        ax.scatter(x[~cross], flin[~cross], s=18, alpha=0.6, color="#36b",
                   label="same side")
        ax.scatter(x[cross], flin[cross], s=18, alpha=0.6, color="#d62",
                   label="cross-barrier")
        sl, ic, _ = linregress(x, flin)
        xs = np.linspace(x.min(), x.max(), 50)
        ax.plot(xs, sl * xs + ic, "k-", lw=1.3)
        ax.set_xlabel(lab)
        ax.set_ylabel(r"$F_{ST}/(1-F_{ST})$")
        ax.set_title(f"genetic vs {lab}\nR$^2$={r2:.2f}")
        ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "exp_c_resistance.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("  wrote figures/exp_c_resistance.png")

    return dict(r2_euclidean=r2_eucl, r2_resistance=r2_res,
                n_pairs=len(tab), n_cross=int(cross.sum()))


# --------------------------------------------------------------------------- #
def main(which=("a", "b", "c")):
    np.random.seed(RNG_SEED)
    t0 = time.time()
    print(f"spaceprime {getattr(sp, '__version__', '?')}  seed={RNG_SEED}  exps={which}")

    # merge into any existing results so individual experiments can be re-run
    path = HERE / "results.json"
    if path.exists():
        results = json.loads(path.read_text())
    else:
        results = {}
    results.setdefault("meta", {}).update(
        dict(spaceprime_version=getattr(sp, "__version__", "?"), seed=RNG_SEED,
             mut_rate=MUT_RATE, recomb_rate=RECOMB_RATE))

    if "a" in which:
        results["experiment_a"] = experiment_a()
    if "b" in which:
        results["experiment_b"] = experiment_b()
    if "c" in which:
        results["experiment_c"] = experiment_c()
    results["meta"]["runtime_sec"] = round(time.time() - t0, 1)

    path.write_text(json.dumps(results, indent=2, default=float))
    print(f"\nTotal runtime {results['meta']['runtime_sec']} s")
    print("wrote results.json")


if __name__ == "__main__":
    import sys

    sel = tuple(sys.argv[1:]) if len(sys.argv) > 1 else ("a", "b", "c")
    main(sel)
