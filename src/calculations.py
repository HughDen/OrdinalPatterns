from itertools import permutations
from collections import defaultdict
from alcoves import region_alcove_generator, address_to_ordinal_pattern
from alcoves import BFS_weak_order_addresses
from math import factorial
from fractions import Fraction

def ordinal_distribution_laplace(n):
    """ 
    Return the probability distribution over all permutation in S_n
    for ordinal patterns in a random walk with steps from a Laplace 
    distribution.

    Keyword arguments:
    n -- Size of the permutations in the returned distribution

        >>> ordinal_distribution_laplace(2)
        {(1, 2): 0.5, (2, 1): 0.5}

    """

    return {p:ordinal_probability_laplace(p) 
            for p in permutations(range(1,n+1))}

def ordinal_distribution_uniform(n):
    """
    Return the probability distribution over all permutations in S_n
    for ordinal patterns in a random walk with steps from a Laplace
    distribution.

    Keyword arguments:
    n -- Size of the permutations in the returned distribution

        >>> ordinal_distribution_uniform(2)
        {(1, 2): 0.5, (2, 1): 0.5}

    """
    distribution = defaultdict(int)
    
    for alcove in region_alcove_generator(n - 1, -1, 0):
        distribution[address_to_ordinal_pattern(alcove)] += 1
    
    return {p:distribution[p] / (2**(n-1) * factorial(n-1))
            for p in permutations(range(1,n+1))}

def ordinal_distribution_uniform_BFS(n):
    """
    Return the probability distribution over all permutations in S_n
    for ordinal patterns in a random walk with steps from a Laplace
    distribution using the breadth first search method
    BFS_weak_order_addresses from alcoves.py.

    Keyword arguments:
    n -- Size of the permutations in the returned distribution

        >>> ordinal_distribution_uniform(2)
        {(1, 2): 0.5, (2, 1): 0.5}

    """
    distribution = {}
    
    for p in permutations(range(1,n+1)):
        distribution[p] = len(BFS_weak_order_addresses(p))
    
    return {p:distribution[p] / (2**(n-1) * factorial(n-1))
            for p in permutations(range(1,n+1))}

def almost_consecutive_bernoulli(n):
    """
    Return the probability that a random walk with n - 1 steps from a
    symmetric density function produces an almost consecutive ordinal
    pattern. Almost consecutive permutations are the ones whose 
    consecutive values are at most two positions apart.

    Keyword arguments:
    n -- The n of S_n

        >>> almost_consecutive_bernoulli(4)
        0.6666666666666666
    """

    bernoulli = 0.0
    for p in permutations(range(1,n+1)):
        if is_almost_consecutive(p):
            bernoulli += ordinal_probability_laplace(p)

    return bernoulli

def almost_consecutive_bernoulli_fraction(n):
    """
    Return the probability that a random walk with n - 1 steps from a
    symmetric density function produces an almost consecutive ordinal
    pattern. Almost consecutive permutations are the ones whose 
    consecutive values are at most two positions apart.

    Keyword arguments:
    n -- The n of S_n

        >>> almost_consecutive_bernoulli_fraction(4)
        Fraction(2, 3)
    """

    bernoulli = Fraction(0)
    for p in permutations(range(1,n+1)):
        if is_almost_consecutive(p):
            bernoulli += _ordinal_probability_laplace_fraction(p)

    return bernoulli

def ordinal_probability_laplace(permutation):
    """ 
    Return the probability that the input permutation occurs as an ordinal
    pattern in a random walk with steps from a Laplace distribution.

    Keyword arguments:
    permutation -- A tuple or list that is a permutation of {1,...,n}

        >>> ordinal_probability_laplace((2,5,1,3,4))
        0.001736111111111111

    """

    N = len(permutation) - 1
    between = between_measure(permutation)
    product = 1

    # As long as the user supplies a permutation, all values in the 
    # between measure tuple are >= 1. Thus, there's no risk of division by zero.

    for i in range(len(between)):
        product = product * (1.0 / between[i])

    return product * (1.0 / 2**N)

def between_measure(permutation):
    """ 
    Return a tuple capturing how often a value is between two values in
    a permutation. Specifically, how often does j in {1,...,n-1} satisfy 
    p(i) <= j < p(i+1) or p(i+1) <= j < p(i)?

    This measure allows us to compute the probability that a random walk
    with steps drawn from the Laplace distribution has the given permutation
    as an ordinal pattern. It also is good for comparing two permutations in
    a similar walk but with Gaussian steps.

    Keyword arguments:
    permutation -- A tuple or list that is a permutation of {1,...,n}

        >>> between_measure((2,5,1,3,4))
        (2, 3, 3, 2)

    """

    between_count = []

    for j in range(len(permutation) - 1):
        count = 0

        # Note that the usage of j + 1 below is based on the assumption that  
        # the user feeds the function a permutation of {1,...,n}

        for i in range(len(permutation) - 1):
            if permutation[i] <= j + 1 < permutation[i + 1] or\
               permutation[i + 1] <= j + 1 < permutation[i]:
                count += 1
        between_count.append(count)

    return tuple(between_count)

def is_almost_consecutive(permutation):
    """
    Returns True if the consecutive values are at most two positions apart.
    Otherwise, it returns False.

    Keyword arguments:
    permutation -- a permutation

        >>> is_almost_consecutive((5,6,4,2,3,1))
        True

    """

    inverse = {permutation[i]:(i+1) for i in range(len(permutation))}
    almost_consecutive = True
    for j in range(2,len(permutation) + 1):
        if abs(inverse[j] - inverse[j-1]) > 2:
            almost_consecutive = False

    return almost_consecutive

def _ordinal_probability_laplace_fraction(permutation):
    """ 
    Return the probability that the input permutation occurs as an ordinal
    pattern in a random walk with steps from a Laplace distribution as a 
    Python Fraction type.

    Keyword arguments:
    permutation -- A tuple or list that is a permutation of {1,...,n}

        >>> _ordinal_probability_laplace_fraction((2,5,1,3,4))
        Fraction(1, 576)

    """

    N = len(permutation) - 1
    between = between_measure(permutation)
    product = Fraction(1)

    # As long as the user supplies a permutation, all values in the 
    # between measure tuple are >= 1. Thus, there's no risk of division by zero.

    for i in range(len(between)):
        product = product * Fraction(1, between[i])

    return product * Fraction(1,2**N)
