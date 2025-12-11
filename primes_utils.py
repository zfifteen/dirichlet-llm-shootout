import math


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
    sieve = [True] * (limit + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            for j in range(i * i, limit + 1, i):
                sieve[j] = False
    return [i for i in range(2, limit + 1) if sieve[i]]


def primes_by_residue_class(limit, modulus):
    """
    Partition primes by residue class modulo N.

    Args:
        limit (int): Upper bound for prime search
        modulus (int): The modulus N for residue classes

    Returns:
        dict: Maps residue class r to list of primes p where p ≡ r (mod N)
    """
    primes = sieve_of_eratosthenes(limit)
    residues = {r: [] for r in range(modulus)}
    for p in primes:
        residues[p % modulus].append(p)
    return residues


def compute_empirical_density(limit, modulus):
    """
    Compute observed frequency of primes in each residue class.

    Args:
        limit (int): Upper bound for analysis
        modulus (int): The modulus N

    Returns:
        dict: Maps residue class r to empirical density (observed frequency)
    """
    residues = primes_by_residue_class(limit, modulus)
    all_primes = sieve_of_eratosthenes(limit)
    total_primes = len(all_primes)
    densities = {}
    coprimes = [r for r in range(modulus) if math.gcd(r, modulus) == 1]
    for r in coprimes:
        count = len(residues[r])
        if total_primes > 0:
            densities[r] = count / total_primes
        else:
            densities[r] = 0
    return densities
