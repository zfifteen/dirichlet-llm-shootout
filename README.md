# Dirichlet's Theorem: Computational Illustration

## Introduction
Dirichlet's theorem guarantees that every arithmetic progression \(a, a+q, a+2q,\dots\) with \(\gcd(a,q)=1\) contains infinitely many primes. Here we build a lightweight numerical exploration to watch those primes distribute themselves across residue classes modulo a chosen \(N\). The goal is to make the phenomenon visible with simple, reproducible code and plots.

Beyond demonstrating convergence to uniform densities among coprime classes, the experiment highlights finite-sample quirks such as prime races where one residue class temporarily leads another. The target audience is a senior undergraduate or early graduate student comfortable with basic number theory and Python.

## Formal Statement
Let \(a\) and \(q\) be positive integers with \(\gcd(a,q)=1\). Dirichlet's theorem states that the arithmetic progression \(a + nq\) contains infinitely many primes. Moreover, primes become equidistributed among the \(\varphi(q)\) residue classes coprime to \(q\): for each such \(a\),
\[
\lim_{x\to\infty} \frac{\pi(x; q, a)}{\pi(x)} = \frac{1}{\varphi(q)},
\]
where \(\pi(x; q, a)\) counts primes \(p\le x\) with \(p\equiv a \pmod q\), and \(\pi(x)\) counts all primes. Our empirical densities approximate these limits for finite \(x\).

## Euler's Totient Function \(\varphi(N)\)
The totient \(\varphi(N)\) is the number of integers in \(\{1,\dots,N\}\) that are coprime to \(N\). It factors as
\[
\varphi(N) = N \prod_{p \mid N} \left(1 - \frac{1}{p}\right),
\]
with the product over distinct prime divisors of \(N\). In Dirichlet's theorem, \(1/\varphi(N)\) is the theoretical density of primes in each coprime residue class modulo \(N\).

## Computational Methodology
1. **Prime generation:** `sieve_of_eratosthenes(limit)` builds all primes up to a chosen cutoff (efficient to at least 1,000,000).
2. **Residue partition:** `primes_by_residue_class(limit, modulus)` bins primes by congruence class modulo \(N\).
3. **Empirical density:** `compute_empirical_density(limit, modulus)` counts primes coprime to \(N\) and reports the share in each residue class (non-coprime residues receive density 0).
4. **Visualization:** `visualizations.py` samples densities at increasing cutoffs, compares empirical and theoretical values, and depicts broader patterns (heatmap and a prime-race plot).

All computations are self-contained—no external data sources are used.

## How to Run
```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\\Scripts\\activate
pip install -r requirements.txt
python visualizations.py
```
PNG files will be written to the repository root: `convergence_N3.png`, `convergence_N4.png`, `convergence_N5.png`, `comparison_N12.png`, `heatmap_multiple_N.png`, and `creative_visualization.png`.

## Interpreting the Results
- **Convergence plots (N=3,4,5):** Empirical densities for each coprime residue drift toward the horizontal line at \(1/\varphi(N)\). Fluctuations shrink as the prime cutoff grows, illustrating Dirichlet's predicted equidistribution.
- **Bar comparison (N=12):** Side-by-side bars show empirical densities at a fixed large cutoff (~200k) against the theoretical \(1/\varphi(12)=1/4\). Small deviations reflect finite-sample noise.
- **Heatmap:** Rows (different \(N\)) and columns (residue classes) display densities; non-coprime classes appear muted. Uniform coloring across coprime columns signals the expected balance.
- **Creative visualization:** A "prime race" for \(N=4\) tracks the count difference between residues 1 and 3. The curve crosses zero repeatedly—evidence that while long-run densities equalize, short-run lead changes are common.

These visuals collectively demonstrate the core phenomenon: primes favor no coprime residue class in the long run, even though local irregularities create temporary biases.
