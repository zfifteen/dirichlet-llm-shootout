# Dirichlet's Theorem on Primes in Arithmetic Progressions

## Introduction

Dirichlet's theorem on arithmetic progressions is a fundamental result in number theory, stating that for any two positive coprime integers \(a\) and \(d\), there are infinitely many prime numbers of the form \(n = a + kd\) where \(k\) is a non-negative integer. This theorem, proved by Peter Gustav Lejeune Dirichlet in 1837, extends the notion of primes being "randomly" distributed while respecting modular constraints.

The theorem implies not only infinitude but also a natural density among all primes. Specifically, the proportion of primes congruent to \(a \pmod{d}\) is \(\frac{1}{\phi(d)}\), where \(\phi\) is Euler's totient function. This computational illustration demonstrates the theorem empirically by generating primes and analyzing their distribution across residue classes modulo various \(N\).

## Formal Statement

**Theorem**: Let \(a\) and \(d > 0\) be coprime integers (i.e., \(\gcd(a, d) = 1\)). Then there are infinitely many primes \(p\) such that \(p \equiv a \pmod{d}\).

Moreover, the natural density of such primes among all primes is \(\frac{1}{\phi(d)}\), where \(\phi(d) = |\{ k : 1 \leq k \leq d, \gcd(k, d) = 1 \}|\) is Euler's totient function, counting the number of integers up to \(d\) that are coprime to \(d\).

## Euler's Totient Function \(\phi(N)\)

The totient function \(\phi(N)\) measures the size of the multiplicative group modulo \(N\). For example:
- \(\phi(3) = 2\) (1, 2 coprime to 3)
- \(\phi(4) = 2\) (1, 3)
- \(\phi(5) = 4\) (1,2,3,4)
- \(\phi(12) = 4\) (1,5,7,11)

It is multiplicative: if \(\gcd(m,n)=1\), then \(\phi(mn) = \phi(m)\phi(n)\).

## Computational Experiment Methodology

We use the Sieve of Eratosthenes to generate all primes up to a large limit (e.g., 100,000). These primes are then partitioned into residue classes modulo \(N\). For each residue class \(r\) coprime to \(N\), we compute the empirical density as the proportion of primes congruent to \(r \pmod{N}\).

Convergence plots show how these densities approach the theoretical value \(\frac{1}{\phi(N)}\) as the limit increases. A comparison for \(N=12\) and a heatmap across multiple \(N\) provide further insights. All visualizations are generated using Matplotlib.

## Instructions to Run

1. Ensure Python 3.9+ is installed.
2. Install dependencies: `pip install -r requirements.txt`
3. Run the visualization script: `python visualizations.py`
4. View the generated PNG files and refer to this README for interpretation.

## Interpreting the Results

The convergence plots (e.g., `convergence_N3.png`) illustrate how empirical densities stabilize towards \(\frac{1}{\phi(N)}\) with increasing prime limits, providing visual evidence of Dirichlet's theorem. Oscillations decrease as more primes are included, confirming the equidistribution.

For \(N=12\), the bar chart `comparison_N12.png` shows empirical densities close to 0.25 for residues 1,5,7,11, matching theory.

The heatmap `heatmap_multiple_N.png` reveals patterns across moduli, with brighter cells indicating higher prime densities in coprime classes.

The creative visualization `creative_visualization.png` plots primes in a specific progression to highlight their infinitude empirically.

This experiment underscores the theorem's predictive power, even for modest computational scales.