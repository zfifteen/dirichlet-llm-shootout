# Dirichlet's Theorem: Primes in Arithmetic Progressions

A computational illustration of one of the most beautiful results in analytic number theory—Dirichlet's theorem on primes in arithmetic progressions.

## Introduction

When we list the prime numbers—2, 3, 5, 7, 11, 13, 17, 19, 23, 29, ...—a natural question arises: how are these primes distributed among different residue classes? For instance, if we look at primes modulo 4, we find that every prime greater than 2 is either congruent to 1 or 3 (mod 4), since primes congruent to 0 or 2 would be divisible by 2. But are there "more" primes of the form 4k+1 or 4k+3? Does one residue class eventually dominate?

In 1837, Peter Gustav Lejeune Dirichlet proved a remarkable theorem: not only are there infinitely many primes in each "allowable" residue class, but these primes are distributed *equally* among all such classes in a precise asymptotic sense. This result, which required the invention of entirely new mathematical machinery (Dirichlet L-functions and characters), stands as one of the foundational achievements of analytic number theory.

This repository provides a computational exploration of Dirichlet's theorem, allowing you to visualize how the empirical distribution of primes converges to the theoretical predictions as we consider larger and larger sets of primes.

## The Theorem

### Formal Statement

**Dirichlet's Theorem on Primes in Arithmetic Progressions:**

Let $N$ and $r$ be positive integers with $\gcd(r, N) = 1$ (i.e., $r$ and $N$ are coprime). Then the arithmetic progression

$$r, \, r + N, \, r + 2N, \, r + 3N, \, \ldots$$

contains infinitely many prime numbers.

Moreover, if $\pi(x; N, r)$ denotes the number of primes $p \leq x$ with $p \equiv r \pmod{N}$, then

$$\lim_{x \to \infty} \frac{\pi(x; N, r)}{\pi(x)} = \frac{1}{\varphi(N)}$$

where $\pi(x)$ is the total count of primes up to $x$, and $\varphi(N)$ is Euler's totient function.

### Interpretation

This asymptotic formula tells us that primes are **equidistributed** among the $\varphi(N)$ residue classes coprime to $N$. Each such class captures, in the limit, exactly $1/\varphi(N)$ of all primes.

For example:
- **Modulo 3:** The coprime residue classes are {1, 2}, so $\varphi(3) = 2$. Each class contains approximately half of all primes.
- **Modulo 4:** The coprime residue classes are {1, 3}, so $\varphi(4) = 2$. Primes of the form $4k+1$ and $4k+3$ are equally common asymptotically.
- **Modulo 10:** The coprime residue classes are {1, 3, 7, 9}, so $\varphi(10) = 4$. Each class contains approximately 25% of primes.

## Euler's Totient Function

Euler's totient function $\varphi(N)$ counts the number of integers $k$ in the range $1 \leq k \leq N$ that are coprime to $N$:

$$\varphi(N) = \#\{k : 1 \leq k \leq N, \, \gcd(k, N) = 1\}$$

Key properties:
- $\varphi(p) = p - 1$ for any prime $p$
- $\varphi(p^k) = p^{k-1}(p-1)$ for prime powers
- $\varphi(mn) = \varphi(m)\varphi(n)$ when $\gcd(m,n) = 1$ (multiplicativity)

| $N$ | Coprime classes | $\varphi(N)$ | Theoretical density |
|-----|-----------------|--------------|---------------------|
| 3   | {1, 2}          | 2            | 1/2 = 0.500         |
| 4   | {1, 3}          | 2            | 1/2 = 0.500         |
| 5   | {1, 2, 3, 4}    | 4            | 1/4 = 0.250         |
| 6   | {1, 5}          | 2            | 1/2 = 0.500         |
| 10  | {1, 3, 7, 9}    | 4            | 1/4 = 0.250         |
| 12  | {1, 5, 7, 11}   | 4            | 1/4 = 0.250         |

## Computational Experiment

### Methodology

Our experiment empirically verifies Dirichlet's theorem by:

1. **Prime Generation:** Using the Sieve of Eratosthenes to generate all primes up to a specified limit (default: 1,000,000).

2. **Classification:** Partitioning primes by their residue modulo $N$ for various values of $N$.

3. **Density Computation:** For each coprime residue class $r$, computing the empirical density:
   $$\hat{d}(r) = \frac{\#\{p \leq x : p \equiv r \pmod{N}\}}{\pi(x)}$$

4. **Convergence Analysis:** Tracking how $\hat{d}(r)$ approaches $1/\varphi(N)$ as the prime limit increases.

### Code Structure

- `primes_utils.py`: Core functions for prime generation and density computation
- `visualizations.py`: Matplotlib-based visualization generation
- `requirements.txt`: Python dependencies

## Running the Code

### Prerequisites

- Python 3.9 or higher
- pip package manager

### Installation

```bash
# Clone or download this gist
cd dirichlet-theorem-illustration

# Install dependencies
pip install -r requirements.txt
```

### Execution

```bash
# Run the utility module demonstration
python primes_utils.py

# Generate all visualizations
python visualizations.py
```

The visualization script will generate six PNG files in the current directory:
- `convergence_N3.png` — Convergence plot for modulus 3
- `convergence_N4.png` — Convergence plot for modulus 4
- `convergence_N5.png` — Convergence plot for modulus 5
- `comparison_N12.png` — Bar chart comparing empirical vs theoretical for modulus 12
- `heatmap_multiple_N.png` — Heatmap across multiple moduli
- `creative_visualization.png` — Error convergence and spiral visualization

## Understanding the Results

### Convergence Plots

The convergence plots (N=3, 4, 5) show empirical densities for each coprime residue class as a function of the prime limit. Key observations:

- **Initial Fluctuations:** For small prime limits, empirical densities show significant variance. This is expected: with few samples, random fluctuations dominate.

- **Asymptotic Convergence:** As the limit increases, all curves converge toward the horizontal line at $1/\varphi(N)$. This beautifully illustrates Dirichlet's theorem in action.

- **Rate of Convergence:** The error decreases approximately as $O(1/\sqrt{\pi(x)})$, consistent with statistical sampling theory.

### Bar Chart Comparison

The bar chart for N=12 provides a snapshot comparison at a fixed limit (1,000,000 primes). The close agreement between empirical (blue) and theoretical (coral) bars demonstrates that even at finite limits, Dirichlet's equidistribution is remarkably accurate.

### Heatmap Analysis

The heatmap visualization shows the ratio of empirical to theoretical density across multiple moduli. Values near 1.0 (white) indicate excellent agreement. This view reveals:

- Consistent equidistribution across different moduli
- The universality of Dirichlet's theorem regardless of the specific modulus chosen

### Creative Visualization

The dual-panel creative visualization offers two perspectives:

1. **Error Convergence (Log-Log Plot):** Shows how the maximum error from the theoretical density decreases as sample size grows. The parallel slopes suggest similar convergence rates across different moduli.

2. **Prime Spiral:** Inspired by the 3Blue1Brown visualization, this polar plot places each prime $p$ at angle $p$ radians and radius $\sqrt{p}$. Primes are colored by their residue modulo 6, revealing the beautiful structure: all primes greater than 3 fall into exactly two residue classes (1 and 5 mod 6).

## Mathematical Significance

Dirichlet's theorem has profound implications:

1. **Density Rather Than Infinity:** The theorem doesn't just say infinitely many primes exist in each class—it quantifies their relative abundance.

2. **Analytic Methods:** Dirichlet's proof pioneered the use of analysis (L-functions, characters) in number theory, establishing a bridge between continuous and discrete mathematics.

3. **Generalization of Euclid:** Euclid proved there are infinitely many primes. Dirichlet showed this infinity is equitably distributed among all "eligible" residue classes.

## References

- Dirichlet, P.G.L. (1837). "Beweis des Satzes, dass jede unbegrenzte arithmetische Progression, deren erstes Glied und Differenz ganze Zahlen ohne gemeinschaftlichen Factor sind, unendlich viele Primzahlen enthält."
- Apostol, T.M. (1976). *Introduction to Analytic Number Theory*. Springer.
- 3Blue1Brown. "Why do prime numbers make these spirals?" [YouTube](https://youtu.be/EK32jo7i5LQ)

## License

This code is provided for educational purposes. Feel free to use, modify, and distribute.

---

*Generated as part of an LLM comparative experiment on mathematical exposition and coding.*
