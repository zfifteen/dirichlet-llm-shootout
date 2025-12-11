"""
Utility functions for working with primes and residue classes.
The sieve is fast enough for limits on the order of 1,000,000+.
"""

from math import gcd, isqrt
from collections import defaultdict


def sieve_of_eratosthenes(limit):
    """
    Generate all primes up to limit using Sieve of Eratosthenes.

    Args:
        limit (int): Upper bound for prime generation

    Returns:
        list: All primes <= limit
    """
    if limit < 2:
        return []
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"  # 0 and 1 are not primes
    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            step = p
            start = p * p
            sieve[start: limit + 1: step] = b"\x00" * (((limit - start) // step) + 1)
    return [i for i, is_prime in enumerate(sieve) if is_prime]


def primes_by_residue_class(limit, modulus):
    """
    Partition primes by residue class modulo N.

    Args:
        limit (int): Upper bound for prime search
        modulus (int): The modulus N for residue classes

    Returns:
        dict: Maps residue class r to list of primes p where p ≡ r (mod N)
    """
    if modulus < 2:
        raise ValueError("modulus must be >= 2")
    primes = sieve_of_eratosthenes(limit)
    classes = {r: [] for r in range(modulus)}
    for p in primes:
        r = p % modulus
        classes[r].append(p)
    return classes


def compute_empirical_density(limit, modulus):
    """
    Compute observed frequency of primes in each residue class.

    Frequencies are computed among primes that are coprime to the modulus,
    so densities across valid residue classes sum to 1. Residues not coprime
    to the modulus are reported with density 0.0.

    Args:
        limit (int): Upper bound for analysis
        modulus (int): The modulus N

    Returns:
        dict: Maps residue class r to empirical density (observed frequency)
    """
    if modulus < 2:
        raise ValueError("modulus must be >= 2")
    primes = sieve_of_eratosthenes(limit)
    coprime_primes = [p for p in primes if gcd(p, modulus) == 1]
    total = len(coprime_primes)
    densities = {r: 0.0 for r in range(modulus)}
    if total == 0:
        return densities
    counts = defaultdict(int)
    for p in coprime_primes:
        counts[p % modulus] += 1
    for r in densities:
        if gcd(r, modulus) == 1:
            densities[r] = counts[r] / total
    return densities


def euler_totient(n):
    """Compute Euler's totient function phi(n)."""
    if n < 1:
        raise ValueError("n must be >= 1")
    result = n
    p = 2
    temp = n
    while p * p <= temp:
        if temp % p == 0:
            while temp % p == 0:
                temp //= p
            result -= result // p
        p += 1
    if temp > 1:
        result -= result // temp
    return result


__all__ = [
    "sieve_of_eratosthenes",
    "primes_by_residue_class",
    "compute_empirical_density",
    "euler_totient",
]
