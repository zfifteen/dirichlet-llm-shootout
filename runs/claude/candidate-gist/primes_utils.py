"""
Prime number utilities for exploring Dirichlet's theorem on primes in arithmetic progressions.

This module provides efficient implementations for:
- Prime generation via the Sieve of Eratosthenes
- Partitioning primes by residue class
- Computing empirical densities for comparison with theoretical predictions
"""

from math import gcd
from typing import Dict, List


def euler_totient(n: int) -> int:
    """
    Compute Euler's totient function φ(n).

    φ(n) counts the number of integers k in [1, n] that are coprime to n,
    i.e., gcd(k, n) = 1.

    Args:
        n (int): Positive integer

    Returns:
        int: The value φ(n)
    """
    if n < 1:
        raise ValueError("n must be a positive integer")

    result = 0
    for k in range(1, n + 1):
        if gcd(k, n) == 1:
            result += 1
    return result


def coprime_residue_classes(n: int) -> List[int]:
    """
    Return all residue classes coprime to n.

    These are the residue classes r in [1, n-1] where gcd(r, n) = 1.
    By Dirichlet's theorem, primes are equidistributed among these classes
    in the limit.

    Args:
        n (int): The modulus

    Returns:
        list: Residue classes coprime to n
    """
    return [r for r in range(1, n) if gcd(r, n) == 1]


def sieve_of_eratosthenes(limit: int) -> List[int]:
    """
    Generate all primes up to limit using Sieve of Eratosthenes.

    The algorithm works by iteratively marking the multiples of each prime
    starting from 2. Numbers that remain unmarked are prime.

    Time complexity: O(n log log n)
    Space complexity: O(n)

    Args:
        limit (int): Upper bound for prime generation

    Returns:
        list: All primes <= limit
    """
    if limit < 2:
        return []

    # Boolean array: is_prime[i] indicates whether i is prime
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False

    # Only need to check up to sqrt(limit)
    i = 2
    while i * i <= limit:
        if is_prime[i]:
            # Mark all multiples of i as composite
            for j in range(i * i, limit + 1, i):
                is_prime[j] = False
        i += 1

    return [i for i in range(limit + 1) if is_prime[i]]


def primes_by_residue_class(limit: int, modulus: int) -> Dict[int, List[int]]:
    """
    Partition primes by residue class modulo N.

    For each prime p <= limit, we compute p mod modulus and group primes
    accordingly. By Dirichlet's theorem, for residue classes r coprime to N,
    there are infinitely many primes congruent to r mod N.

    Args:
        limit (int): Upper bound for prime search
        modulus (int): The modulus N for residue classes

    Returns:
        dict: Maps residue class r to list of primes p where p ≡ r (mod N)
    """
    if modulus < 1:
        raise ValueError("modulus must be a positive integer")

    primes = sieve_of_eratosthenes(limit)

    # Initialize dictionary for all residue classes
    residue_to_primes: Dict[int, List[int]] = {r: [] for r in range(modulus)}

    for p in primes:
        residue = p % modulus
        residue_to_primes[residue].append(p)

    return residue_to_primes


def compute_empirical_density(limit: int, modulus: int) -> Dict[int, float]:
    """
    Compute observed frequency of primes in each residue class.

    The empirical density for residue class r is:
        (# of primes p ≤ limit with p ≡ r mod N) / (total # of primes ≤ limit)

    By Dirichlet's theorem, for residue classes r coprime to N, this density
    converges to 1/φ(N) as limit → ∞.

    Args:
        limit (int): Upper bound for analysis
        modulus (int): The modulus N

    Returns:
        dict: Maps residue class r to empirical density (observed frequency)
    """
    residue_to_primes = primes_by_residue_class(limit, modulus)

    # Count total primes (exclude 0 residue class if modulus divides a prime)
    total_primes = sum(len(primes) for primes in residue_to_primes.values())

    if total_primes == 0:
        return {r: 0.0 for r in range(modulus)}

    densities: Dict[int, float] = {}
    for residue, primes in residue_to_primes.items():
        densities[residue] = len(primes) / total_primes

    return densities


def theoretical_density(modulus: int) -> float:
    """
    Compute the theoretical asymptotic density for primes in coprime residue classes.

    By Dirichlet's theorem, primes are equidistributed among the φ(N) residue
    classes coprime to N, so each such class has asymptotic density 1/φ(N).

    Args:
        modulus (int): The modulus N

    Returns:
        float: The theoretical density 1/φ(N)
    """
    return 1.0 / euler_totient(modulus)


def prime_counts_up_to_limits(limits: List[int], modulus: int) -> Dict[int, List[int]]:
    """
    Compute cumulative prime counts in each residue class for multiple limits.

    This is useful for visualizing convergence to the theoretical density.

    Args:
        limits (list): List of upper bounds to analyze
        modulus (int): The modulus N

    Returns:
        dict: Maps residue class r to list of cumulative counts at each limit
    """
    max_limit = max(limits)
    primes = sieve_of_eratosthenes(max_limit)

    # Precompute residue for each prime
    prime_residues = [(p, p % modulus) for p in primes]

    # Initialize counts
    coprime_classes = coprime_residue_classes(modulus)
    counts: Dict[int, List[int]] = {r: [] for r in coprime_classes}

    for limit in sorted(limits):
        # Count primes up to this limit in each coprime residue class
        for r in coprime_classes:
            count = sum(1 for p, res in prime_residues if p <= limit and res == r)
            counts[r].append(count)

    return counts


def densities_at_limits(limits: List[int], modulus: int) -> Dict[int, List[float]]:
    """
    Compute empirical densities at multiple limits for convergence analysis.

    Args:
        limits (list): List of upper bounds to analyze
        modulus (int): The modulus N

    Returns:
        dict: Maps coprime residue class r to list of densities at each limit
    """
    max_limit = max(limits)
    primes = sieve_of_eratosthenes(max_limit)

    coprime_classes = coprime_residue_classes(modulus)
    densities: Dict[int, List[float]] = {r: [] for r in coprime_classes}

    for limit in sorted(limits):
        primes_up_to_limit = [p for p in primes if p <= limit]
        total = len(primes_up_to_limit)

        if total == 0:
            for r in coprime_classes:
                densities[r].append(0.0)
            continue

        for r in coprime_classes:
            count = sum(1 for p in primes_up_to_limit if p % modulus == r)
            densities[r].append(count / total)

    return densities


if __name__ == "__main__":
    # Quick demonstration
    print("Dirichlet's Theorem: Prime Distribution in Residue Classes")
    print("=" * 60)

    limit = 100000

    for N in [3, 4, 5, 10, 12]:
        phi_N = euler_totient(N)
        theoretical = theoretical_density(N)
        empirical = compute_empirical_density(limit, N)
        coprime = coprime_residue_classes(N)

        print(f"\nModulus N = {N}, φ(N) = {phi_N}")
        print(f"Theoretical density: 1/φ(N) = {theoretical:.6f}")
        print("Empirical densities for coprime residue classes:")

        for r in coprime:
            print(f"  r = {r}: {empirical[r]:.6f}")