"""
Visualizations for Dirichlet's Theorem on Primes in Arithmetic Progressions.

This module generates six publication-quality figures illustrating how primes
distribute themselves among residue classes, demonstrating convergence to
the theoretical predictions of Dirichlet's theorem.

Output files:
    - convergence_N3.png
    - convergence_N4.png
    - convergence_N5.png
    - comparison_N12.png
    - heatmap_multiple_N.png
    - creative_visualization.png
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from typing import List, Dict

from primes_utils import (
    sieve_of_eratosthenes,
    euler_totient,
    coprime_residue_classes,
    theoretical_density,
    densities_at_limits,
    compute_empirical_density,
    primes_by_residue_class,
)

# Style configuration for publication-quality plots
plt.rcParams.update({
    'font.size': 11,
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'legend.fontsize': 10,
    'figure.dpi': 150,
    'savefig.dpi': 150,
    'savefig.bbox': 'tight',
})


def generate_convergence_limits(max_limit: int = 1_000_000, num_points: int = 50) -> List[int]:
    """Generate logarithmically spaced limits for convergence analysis."""
    return [int(x) for x in np.logspace(3, np.log10(max_limit), num_points)]


def plot_convergence(modulus: int, filename: str, max_limit: int = 1_000_000):
    """
    Generate a convergence plot showing how empirical densities approach
    the theoretical value 1/φ(N) as the prime limit increases.

    Args:
        modulus: The modulus N for residue classes
        filename: Output filename for the PNG
        max_limit: Maximum prime limit to consider
    """
    limits = generate_convergence_limits(max_limit)
    densities = densities_at_limits(limits, modulus)
    theoretical = theoretical_density(modulus)
    phi_N = euler_totient(modulus)
    coprime = coprime_residue_classes(modulus)

    fig, ax = plt.subplots(figsize=(10, 6))

    # Color palette for residue classes
    colors = plt.cm.tab10(np.linspace(0, 1, len(coprime)))

    for i, r in enumerate(coprime):
        ax.plot(limits, densities[r], label=f'r ≡ {r} (mod {modulus})',
                color=colors[i], linewidth=1.5, alpha=0.8)

    # Theoretical horizontal line
    ax.axhline(y=theoretical, color='red', linestyle='--', linewidth=2,
               label=f'Theoretical: 1/φ({modulus}) = 1/{phi_N} ≈ {theoretical:.4f}')

    ax.set_xscale('log')
    ax.set_xlabel('Upper Limit of Primes Considered')
    ax.set_ylabel('Empirical Density (Fraction of Primes)')
    ax.set_title(f'Convergence of Prime Density in Residue Classes mod {modulus}\n'
                 f'(Dirichlet\'s Theorem: Expected density = 1/φ({modulus}) = 1/{phi_N})')
    ax.legend(loc='best', framealpha=0.9)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, max(0.6, theoretical * 2))

    plt.tight_layout()
    plt.savefig(filename)
    plt.close()
    print(f"Saved: {filename}")


def plot_comparison_bar_chart(modulus: int, filename: str, limit: int = 1_000_000):
    """
    Generate a bar chart comparing empirical vs theoretical densities
    for all coprime residue classes.

    Args:
        modulus: The modulus N
        filename: Output filename
        limit: Prime limit for empirical computation
    """
    empirical = compute_empirical_density(limit, modulus)
    theoretical = theoretical_density(modulus)
    phi_N = euler_totient(modulus)
    coprime = coprime_residue_classes(modulus)

    fig, ax = plt.subplots(figsize=(10, 6))

    x = np.arange(len(coprime))
    width = 0.35

    empirical_vals = [empirical[r] for r in coprime]
    theoretical_vals = [theoretical] * len(coprime)

    bars1 = ax.bar(x - width/2, empirical_vals, width, label='Empirical',
                   color='steelblue', alpha=0.8)
    bars2 = ax.bar(x + width/2, theoretical_vals, width, label='Theoretical',
                   color='coral', alpha=0.8)

    ax.set_xlabel('Residue Class r (coprime to N)')
    ax.set_ylabel('Density (Fraction of Primes)')
    ax.set_title(f'Prime Distribution in Residue Classes mod {modulus}\n'
                 f'Primes up to {limit:,} | φ({modulus}) = {phi_N} | '
                 f'Theoretical = 1/{phi_N} ≈ {theoretical:.4f}')
    ax.set_xticks(x)
    ax.set_xticklabels([f'r = {r}' for r in coprime])
    ax.legend(loc='upper right')
    ax.grid(True, axis='y', alpha=0.3)

    # Add value labels on bars
    for bar, val in zip(bars1, empirical_vals):
        height = bar.get_height()
        ax.annotate(f'{val:.4f}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3), textcoords="offset points",
                    ha='center', va='bottom', fontsize=8)

    plt.tight_layout()
    plt.savefig(filename)
    plt.close()
    print(f"Saved: {filename}")


def plot_heatmap_multiple_N(filename: str, moduli: List[int] = None, limit: int = 1_000_000):
    """
    Generate a heatmap showing prime counts across different moduli and
    their coprime residue classes.

    Args:
        filename: Output filename
        moduli: List of moduli to analyze
        limit: Prime limit
    """
    if moduli is None:
        moduli = [3, 4, 5, 6, 8, 10, 12]

    # Find maximum number of coprime residue classes
    max_coprime = max(euler_totient(N) for N in moduli)

    # Build matrix: rows = moduli, columns = residue classes (padded)
    # Value = normalized density (empirical / theoretical)
    matrix = np.full((len(moduli), max_coprime), np.nan)
    row_labels = []
    col_positions = []

    for i, N in enumerate(moduli):
        empirical = compute_empirical_density(limit, N)
        theoretical = theoretical_density(N)
        coprime = coprime_residue_classes(N)
        phi_N = euler_totient(N)

        row_labels.append(f'N={N} (φ={phi_N})')

        for j, r in enumerate(coprime):
            # Store ratio of empirical to theoretical (1.0 = perfect agreement)
            if theoretical > 0:
                matrix[i, j] = empirical[r] / theoretical
            col_positions.append((i, j, r))

    fig, ax = plt.subplots(figsize=(12, 7))

    # Custom colormap: blue (below theoretical) - white (at theoretical) - red (above)
    cmap = plt.cm.RdYlBu_r

    im = ax.imshow(matrix, cmap=cmap, aspect='auto', vmin=0.8, vmax=1.2)

    # Colorbar
    cbar = ax.figure.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Empirical / Theoretical Density Ratio\n(1.0 = Perfect Agreement)', fontsize=10)

    # Labels
    ax.set_yticks(np.arange(len(moduli)))
    ax.set_yticklabels(row_labels)
    ax.set_xlabel('Coprime Residue Class Index')
    ax.set_ylabel('Modulus N')
    ax.set_title(f'Prime Distribution Heatmap: Empirical vs Theoretical Density\n'
                 f'(Primes up to {limit:,})')

    # Add text annotations for valid cells
    for i, N in enumerate(moduli):
        coprime = coprime_residue_classes(N)
        for j, r in enumerate(coprime):
            if not np.isnan(matrix[i, j]):
                text_color = 'white' if abs(matrix[i, j] - 1.0) > 0.1 else 'black'
                ax.text(j, i, f'{r}\n{matrix[i, j]:.2f}',
                        ha='center', va='center', color=text_color, fontsize=8)

    plt.tight_layout()
    plt.savefig(filename)
    plt.close()
    print(f"Saved: {filename}")


def plot_creative_visualization(filename: str, limit: int = 100_000):
    """
    Create a polar/spiral visualization showing primes colored by residue class.

    This visualization is inspired by the 3Blue1Brown video, showing how
    primes form patterns when plotted on a spiral with different moduli.

    Args:
        filename: Output filename
        limit: Maximum prime to include
    """
    primes = sieve_of_eratosthenes(limit)
    modulus = 6  # Creates a nice pattern: primes > 3 are ≡ 1 or 5 (mod 6)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left plot: Error convergence for multiple moduli
    ax1 = axes[0]
    test_moduli = [3, 4, 5, 6, 7, 8, 10, 12]
    limits = generate_convergence_limits(limit, 40)

    for N in test_moduli:
        densities = densities_at_limits(limits, N)
        theoretical = theoretical_density(N)
        coprime = coprime_residue_classes(N)

        # Compute max absolute error across all coprime classes
        max_errors = []
        for lim_densities in zip(*[densities[r] for r in coprime]):
            max_error = max(abs(d - theoretical) for d in lim_densities)
            max_errors.append(max_error)

        ax1.plot(limits, max_errors, label=f'N={N}', linewidth=1.5, alpha=0.8)

    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('Upper Limit of Primes')
    ax1.set_ylabel('Maximum Absolute Error from 1/φ(N)')
    ax1.set_title('Convergence Rate: Error Decreases as Sample Size Grows')
    ax1.legend(loc='best', ncol=2, fontsize=9)
    ax1.grid(True, alpha=0.3)

    # Right plot: Spiral visualization
    ax2 = axes[1]

    # Use primes up to a smaller limit for visual clarity
    spiral_limit = 10000
    spiral_primes = [p for p in primes if p <= spiral_limit]

    # Spiral coordinates: prime p is at angle p and radius sqrt(p)
    angles = np.array(spiral_primes)
    radii = np.sqrt(np.array(spiral_primes))

    # Color by residue mod 6
    colors_map = {1: 'blue', 5: 'red'}
    colors = [colors_map.get(p % modulus, 'gray') for p in spiral_primes]

    ax2.scatter(angles, radii, c=colors, s=3, alpha=0.6)
    ax2.set_xlabel('Angle (radians, equal to prime value)')
    ax2.set_ylabel('Radius (√p)')
    ax2.set_title(f'Prime Spiral (primes ≤ {spiral_limit})\n'
                  f'Blue: p ≡ 1 (mod 6), Red: p ≡ 5 (mod 6)')

    # Convert to polar-like appearance by using actual polar projection
    # We'll create a proper polar subplot instead
    ax2.remove()
    ax2 = fig.add_subplot(122, projection='polar')

    # For polar plot, angle should be in radians, scaled appropriately
    # Use p radians (mod 2π) for angle, sqrt(p) for radius
    angles_polar = np.array(spiral_primes) % (2 * np.pi * 10)  # Multiple rotations
    radii_polar = np.sqrt(np.array(spiral_primes))

    scatter = ax2.scatter(spiral_primes, radii_polar, c=colors, s=2, alpha=0.5)
    ax2.set_title(f'Prime Spiral Visualization\n'
                  f'Angle = p, Radius = √p', pad=20)
    ax2.set_rticks([])  # Hide radial ticks for cleaner look

    # Add legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='blue',
               markersize=8, label='p ≡ 1 (mod 6)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='red',
               markersize=8, label='p ≡ 5 (mod 6)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
               markersize=8, label='p = 2, 3'),
    ]
    ax2.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.3, 1.0))

    plt.tight_layout()
    plt.savefig(filename)
    plt.close()
    print(f"Saved: {filename}")


def generate_all_visualizations():
    """Generate all required visualizations for the Dirichlet theorem illustration."""
    print("Generating visualizations for Dirichlet's Theorem...")
    print("=" * 60)

    # Convergence plots for N = 3, 4, 5
    plot_convergence(3, 'convergence_N3.png')
    plot_convergence(4, 'convergence_N4.png')
    plot_convergence(5, 'convergence_N5.png')

    # Comparison bar chart for N = 12
    plot_comparison_bar_chart(12, 'comparison_N12.png')

    # Heatmap for multiple N values
    plot_heatmap_multiple_N('heatmap_multiple_N.png')

    # Creative visualization
    plot_creative_visualization('creative_visualization.png')

    print("=" * 60)
    print("All visualizations generated successfully!")


if __name__ == "__main__":
    generate_all_visualizations()
