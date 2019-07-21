from math import log, sqrt, exp
from itertools import permutations
from scipy.stats import entropy

def chebyshev_bound(outcome, mean, sd):
    """
    Return the maximum probability guaranteed by Chebyshev's inequality
    that |outcome - mean| is as high (or low) as it is. Chebyshev's inequality
    applies to any distribution.

    Keyword arguments:
    outcome -- Outcome of an experiment
    mean -- mean of the hypothesized distribution
    sd -- standard deviation of the hypothesized distribution

        This example calculates the bound in the case of 75 out of 100 fair
        coin tosses turning up heads.

        >>> chebyshev_bound(65,50, 5)
        0.1111111111111111

    """

    sd_distance = abs(outcome - mean) / sd
    return 1 / sd_distance**2

def hoeffding_chernoff_bound(successes, n, p):
    """
    Return the maximum probability of the outcome of a Bernoulli trial
    guaranteed by a bound due to Hoeffding based on Chernoff bounds. 

    Keyword arguments:
    successes -- Number of successes (or heads or whatever) in a Bernoulli trial
    n -- number of trials
    p -- probability of success

        This example calculates the bound in the case of 75 out of 100 fair
        coin tosses turning up heads.

        >>> hoeffding_chernoff_bound(65,100, 0.5)
        0.02071479764612625

    """

    epsilon = (successes / n) - p
    p_distribution = {'H' : p, 'T' : 1 - p}
    plus_distribution = {'H' : p + epsilon, 'T' : 1 - p - epsilon}
    minus_distribution = {'H' : p - epsilon, 'T' : 1 - p + epsilon}

    kl_plus = kullback_leibler(plus_distribution, p_distribution)
    kl_minus = kullback_leibler(minus_distribution, p_distribution)

    return exp(-kl_plus * n) + exp(-kl_minus * n)

def total_variation_distance(distribution1, distribution2):
    """
    Return the total variation between two discrete distributions. This
    is defined as 1/2 of the sum of the absolute difference between the values.
    The constant 1/2 multiplier enforces a value between 0 and 1 for the 
    intended application of probability distributions.

    Keyword arguments:
    distribution1 -- A dictionary that maps to a numeric type  
    distribution2 -- A dictionary that maps to a numeric type

    A key error will occur if the distributions don't have the same set of keys.

        >>> A = {(1,2):0.3, (2,1):0.7}
        >>> B = {(1,2):0.1, (2,1):0.9}
        >>> total_variation_distance(A,B)
        0.2

    """

    return sum([0.5*abs(distribution1[k] - distribution2[k])
                for k in distribution1.keys()])

def kullback_leibler(distribution1, distribution2):
    """
    Return the kullback-leibler divergence between two distributions. 

    Since scipy.stats.entropy already calculates kl-divergence, this
    just converts the dictionaries to aligned lists before calling entropy.

    Keyword arguments:
    distribution1 -- A dictionary that maps to a numeric type  
    distribution2 -- A dictionary that maps to a numeric type

    A key error can occur if the distributions don't have the same set of keys.

        >>> A = {(1,2):0.3, (2,1):0.7}
        >>> B = {(1,2):0.1, (2,1):0.9}
        >>> kullback_leibler(A,B)
        0.15366358680379857

    """

    # The iteration order of the keys in the two dictionaries need not match.
    # Thus, we align them in their respective lists before calling entropy.

    list1 = []
    list2 = []
    for k in distribution1.keys():
        list1.append(distribution1[k])
        list2.append(distribution2[k])

    return entropy(list1,list2)

def create_bernoulli(distribution1, distribution2):
    """
    Return a tuple of permutations and probabilities p,q such that the 
    permutations have probability p of occurring in distribution1 and
    probability q of occurring in distribution2. The result can be viewed
    as a Bernoulli trial for whether an ordinal pattern of the tuple occurs
    for each distribution. 

    Keyword arguments:
    distribution1 -- A dictionary that maps permutations to a numeric type  
    distribution2 -- A dictionary that maps permutations to a numeric type

    A key error can occur if the distributions don't have the same set of keys.

        >>> A = {(1,2):0.7, (2,1):0.3}
        >>> B = {(1,2):0.4, (2,1):0.6}
        >>> create_bernoulli(A,B)
        (((1, 2),), 0.7, 0.4)

    """

    permutations_list = []

    differences = [(distribution1[k] - distribution2[k], k)
                   for k in distribution1.keys()]
    
    sorted_differences = sorted(differences, reverse=True)

    # The idea is to maximize the difference between the two distributions
    # while keeping one of the created Bernoulli trials near 0.5.
    # This algorithm may not be optimal, but it is straightforward: collect
    # the differences, sort them from highest to lowest, and then collect
    # the associated permutations stopping after the total sum of one
    # distribution exceeds 0.5.
    
    tally = 0.0
    threshold = 0.5
    
    for d,p in sorted_differences:
        tally += distribution1[p]
        permutations_list.append(p)
        if tally > 0.5:
            break

    return tuple(permutations_list),\
           sum([distribution1[p] for p in permutations_list]),\
           sum([distribution2[p] for p in permutations_list])
