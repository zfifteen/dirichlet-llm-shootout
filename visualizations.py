import matplotlib.pyplot as plt
import numpy as np
from primes_utils import sieve_of_eratosthenes, compute_empirical_density
import math

MAX_LIMIT = 100000
STEP = 5000


def plot_convergence(N, max_limit=MAX_LIMIT, step=STEP):
    primes = sieve_of_eratosthenes(max_limit)
    coprimes = [r for r in range(N) if math.gcd(r, N) == 1]
    phi = len(coprimes)
    theoretical = 1 / phi

    limits = list(range(step, max_limit + 1, step))
    densities_dict = {r: [] for r in coprimes}

    cumulative_primes = []
    for lim in limits:
        p_up_to_lim = [p for p in primes if p <= lim]
        cumulative_primes.append(p_up_to_lim)

    for r in coprimes:
        for p_list in cumulative_primes:
            total = len(p_list)
            count = sum(1 for p in p_list if p % N == r)
            densities_dict[r].append(count / total if total > 0 else 0)

    fig, ax = plt.subplots(figsize=(10, 6))
    for r in coprimes:
        ax.plot(limits, densities_dict[r], label=f"r={r} mod {N}")
    ax.axhline(
        y=theoretical,
        color="k",
        linestyle="--",
        label=f"Theoretical: 1/φ({N}) = {theoretical}",
    )
    ax.set_xlabel("Upper Limit of Primes")
    ax.set_ylabel("Empirical Density")
    ax.set_title(f"Convergence of Prime Densities for Modulus N={N}")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"convergence_N{N}.png", dpi=300)
    plt.close()


# Generate convergence plots
plot_convergence(3)
plot_convergence(4)
plot_convergence(5)


# Comparison for N=12
def plot_comparison_N12(limit=MAX_LIMIT):
    coprimes_12 = [1, 5, 7, 11]
    phi_12 = 4
    theoretical = 1 / phi_12
    emp_densities = compute_empirical_density(limit, 12)

    r_values = coprimes_12
    emp_values = [emp_densities[r] for r in r_values]
    theo_values = [theoretical] * len(r_values)

    x = np.arange(len(r_values))
    width = 0.35

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(x - width / 2, emp_values, width, label="Empirical", alpha=0.8)
    ax.bar(x + width / 2, theo_values, width, label="Theoretical", alpha=0.8)

    ax.set_xlabel("Residue Class r mod 12")
    ax.set_ylabel("Density")
    ax.set_title("Empirical vs Theoretical Densities for N=12")
    ax.set_xticks(x)
    ax.set_xticklabels([str(r) for r in r_values])
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("comparison_N12.png", dpi=300)
    plt.close()


plot_comparison_N12()


# Heatmap for multiple N
def plot_heatmap_multiple_N(max_limit=MAX_LIMIT):
    Ns = [3, 4, 5, 6, 8, 10, 12]
    max_N = 12
    data = np.zeros((len(Ns), max_N))

    for i, N in enumerate(Ns):
        coprimes = [r for r in range(N) if math.gcd(r, N) == 1]
        densities = compute_empirical_density(max_limit, N)
        for r in coprimes:
            data[i, r] = densities[r] * len(
                coprimes
            )  # Normalize to show relative density

    fig, ax = plt.subplots(figsize=(12, 8))
    im = ax.imshow(data, cmap="YlOrRd", aspect="auto")
    ax.set_xticks(range(max_N))
    ax.set_yticks(range(len(Ns)))
    ax.set_yticklabels([str(n) for n in Ns])
    ax.set_xlabel("Residue Class r")
    ax.set_ylabel("Modulus N")
    ax.set_title("Heatmap of Normalized Prime Densities Across Moduli")
    plt.colorbar(im, ax=ax, label="Normalized Density")
    plt.tight_layout()
    plt.savefig("heatmap_multiple_N.png", dpi=300)
    plt.close()


plot_heatmap_multiple_N()


# Creative visualization: Primes in arithmetic progression for N=5, r=1
def plot_creative_visualization(limit=MAX_LIMIT):
    N = 5
    r = 1
    primes = sieve_of_eratosthenes(limit)
    ap_primes = [p for p in primes if p % N == r and p > N]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    # Plot 1: Primes in AP on log scale
    ax1.scatter(np.log(ap_primes), range(len(ap_primes)), s=10, alpha=0.6)
    ax1.set_xlabel("log(p)")
    ax1.set_ylabel("Index in Sequence")
    ax1.set_title(f"Primes ≡ {r} mod {N} (Log Scale Distribution)")
    ax1.grid(True, alpha=0.3)

    # Plot 2: Gaps between consecutive AP primes
    gaps = np.diff(ap_primes)
    ax2.plot(gaps, marker="o", linestyle="-", markersize=3)
    ax2.set_xlabel("Index")
    ax2.set_ylabel("Gap Size")
    ax2.set_title(f"Gaps in Primes ≡ {r} mod {N}")
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("creative_visualization.png", dpi=300)
    plt.close()


plot_creative_visualization()

print("All visualizations generated successfully.")
