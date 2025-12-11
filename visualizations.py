"""
Generate visualizations illustrating Dirichlet's theorem on primes in arithmetic progressions.
Run `python visualizations.py` to produce PNGs in the current directory.
"""

import math
from collections import defaultdict

import matplotlib

# Use a non-interactive backend so the script works headlessly
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from primes_utils import (
    sieve_of_eratosthenes,
    compute_empirical_density,
    euler_totient,
)


# --------- Helper utilities ---------

def densities_over_limits(modulus, limits):
    """Compute empirical densities at multiple cutoff values.

    Args:
        modulus (int): modulus N
        limits (list[int]): ascending list of cutoffs

    Returns:
        list[dict]: densities per residue for each limit
    """
    max_limit = max(limits)
    primes = sieve_of_eratosthenes(max_limit)
    limits_iter = iter(sorted(limits))
    current_limit = next(limits_iter)
    results = []
    counts = defaultdict(int)
    total_coprime = 0

    def snapshot():
        densities = {r: 0.0 for r in range(modulus)}
        if total_coprime == 0:
            return densities
        for r in range(modulus):
            if math.gcd(r, modulus) == 1:
                densities[r] = counts[r] / total_coprime
        return densities

    for p in primes:
        if math.gcd(p, modulus) == 1:
            total_coprime += 1
            counts[p % modulus] += 1
        while p >= current_limit:
            results.append(snapshot())
            try:
                current_limit = next(limits_iter)
            except StopIteration:
                return results
    # in case last limit is larger than largest prime <= max_limit
    while True:
        results.append(snapshot())
        try:
            current_limit = next(limits_iter)
        except StopIteration:
            break
    return results


def plot_convergence(modulus, filename, limits=None):
    if limits is None:
        limits = [10_000, 20_000, 40_000, 80_000, 120_000, 200_000]
    densities_list = densities_over_limits(modulus, limits)
    phi_val = euler_totient(modulus)
    theoretical = 1 / phi_val

    plt.figure(figsize=(8, 5))
    residue_classes = [r for r in range(modulus) if math.gcd(r, modulus) == 1]
    for r in residue_classes:
        y = [dens[r] for dens in densities_list]
        plt.plot(limits, y, marker="o", label=f"r ≡ {r} mod {modulus}")
    plt.axhline(theoretical, color="k", linestyle="--", linewidth=1.2,
                label=f"theoretical 1/φ({modulus})={theoretical:.3f}")
    plt.title(f"Empirical density convergence (N={modulus})")
    plt.xlabel("Prime upper limit")
    plt.ylabel("Empirical density")
    plt.ylim(bottom=0)
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename, dpi=200)
    plt.close()


def plot_comparison_bar(modulus, limit, filename):
    densities = compute_empirical_density(limit, modulus)
    residue_classes = [r for r in range(modulus) if math.gcd(r, modulus) == 1]
    empirical = [densities[r] for r in residue_classes]
    phi_val = euler_totient(modulus)
    theoretical = [1 / phi_val for _ in residue_classes]

    x = np.arange(len(residue_classes))
    width = 0.35
    plt.figure(figsize=(8, 5))
    plt.bar(x - width / 2, empirical, width, label="Empirical")
    plt.bar(x + width / 2, theoretical, width, label="Theoretical")
    plt.xticks(x, [f"{r} mod {modulus}" for r in residue_classes])
    plt.ylabel("Density")
    plt.title(f"Empirical vs theoretical densities (N={modulus}, limit={limit})")
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename, dpi=200)
    plt.close()


def plot_heatmap(moduli, limit, filename):
    max_residue = max(moduli)
    data = np.full((len(moduli), max_residue), np.nan)
    for i, N in enumerate(moduli):
        densities = compute_empirical_density(limit, N)
        for r in range(N):
            if math.gcd(r, N) == 1:
                data[i, r] = densities[r]
    plt.figure(figsize=(10, 4 + 0.3 * len(moduli)))
    cmap = plt.cm.viridis
    cmap.set_bad(color="lightgray")
    img = plt.imshow(data, aspect="auto", cmap=cmap, interpolation="nearest")
    plt.yticks(range(len(moduli)), [f"N={N}" for N in moduli])
    plt.xticks(range(max_residue), [str(r) for r in range(max_residue)])
    plt.xlabel("Residue class")
    plt.title(f"Prime densities by residue class (limit={limit})")
    cbar = plt.colorbar(img)
    cbar.set_label("Empirical density")
    plt.tight_layout()
    plt.savefig(filename, dpi=200)
    plt.close()


def plot_prime_race(modulus, residues, limit, filename):
    """Creative plot: prime number race between residue classes."""
    primes = sieve_of_eratosthenes(limit)
    race_counts = {r: 0 for r in residues}
    xs, deltas = [], []
    for p in primes:
        if math.gcd(p, modulus) != 1:
            continue
        r = p % modulus
        if r in race_counts:
            race_counts[r] += 1
        delta = race_counts[residues[0]] - race_counts[residues[1]]
        xs.append(p)
        deltas.append(delta)
    plt.figure(figsize=(8, 5))
    plt.plot(xs, deltas, color="#d1495b", linewidth=1.4)
    plt.axhline(0, color="k", linestyle="--", linewidth=1)
    plt.title(f"Prime race: {residues[0]} vs {residues[1]} (mod {modulus})")
    plt.xlabel("Prime upper limit")
    plt.ylabel(f"Count difference (r≡{residues[0]} - r≡{residues[1]})")
    plt.tight_layout()
    plt.savefig(filename, dpi=200)
    plt.close()


# --------- Main script ---------

def main():
    # Convergence plots
    plot_convergence(3, "convergence_N3.png")
    plot_convergence(4, "convergence_N4.png")
    plot_convergence(5, "convergence_N5.png")

    # Comparison bar chart for N=12
    plot_comparison_bar(12, limit=200_000, filename="comparison_N12.png")

    # Heatmap over several moduli
    moduli = [3, 4, 5, 6, 8, 10, 12]
    plot_heatmap(moduli, limit=200_000, filename="heatmap_multiple_N.png")

    # Creative visualization: prime race mod 4 between residues 1 and 3
    plot_prime_race(4, residues=(1, 3), limit=200_000, filename="creative_visualization.png")


if __name__ == "__main__":
    main()
