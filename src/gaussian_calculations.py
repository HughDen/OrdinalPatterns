from math import pi, asin
from itertools import permutations

def gaussian4(H):
    """ 
    Return a dictionary mapping permutations to the probability they occur
    as an ordinal patter in a Brownian motion with Hurst parameter H.
    """

    probabilities = _gaussian4_dictionary(H)

    for p in permutations(range(1,5)):
        if p not in probabilities.keys():
            if _reverse(p) in probabilities.keys():
                probabilities[p] = probabilities[_reverse(p)]
            if _complement(p) in probabilities.keys():
                probabilities[p] = probabilities[_complement(p)]
            if _reverse_complement(p) in probabilities.keys():
                probabilities[p] = probabilities[_reverse_complement(p)]

    return probabilities

def _alpha1(H):
    """
    Return alpha_1 parameter from Bandt and Shiha's Order Patterns in Time 
    Series paper for fractional Brownian motion.
    
    Keyword arguments:
    H -- Hurst parameter H; H = 1/2 is ordinary brownian motion

        >>> _alpha1(0.5)
        0.0

    """

    return (1 + 3**(2*H) - 2**(2*H+1)) / 2

def _alpha2(H):
    """
    Return alpha_2 parameter from Bandt and Shiha's Order Patterns in Time 
    Series paper for fractional Brownian motion.
    
    Keyword arguments:
    H -- Hurst parameter H; H = 1/2 is ordinary brownian motion

        >>> _alpha2(0.5)
        0.0

    """

    return 2**(2*H-1) - 1

def _alpha3(H):
    """
    Return alpha_3 parameter from Bandt and Shiha's Order Patterns in Time Series paper
    for fractional Brownian motion.
    
    Keyword arguments:
    H -- Hurst parameter H; H = 1/2 is ordinary brownian motion

        >>> _alpha3(0.5)
        -0.8164965809277261

    """

    return (1 - 3**(2*H) - 2**(2*H)) / (2*6**H)

def _alpha4(H):
    """
    Return alpha_4 parameter from Bandt and Shiha's Order Patterns in Time 
    Series paper for fractional Brownian motion.
    
    Keyword arguments:
    H -- Hurst parameter H; H = 1/2 is ordinary brownian motion

        >>> _alpha4(0.5)
        0.5

    """

    return (3**(2*H) - 1) / 2**(2*H + 1)

def _alpha5(H):
    """
    Return alpha_5 parameter from Bandt and Shiha's Order Patterns in Time 
    Series paper for fractional Brownian motion.
    
    Keyword arguments:
    H -- Hurst parameter H; H = 1/2 is ordinary brownian motion

        >>> _alpha5(0.5)
        0.7071067811865476

    """

    return 2**(H-1)

def _alpha6(H):
    """
    Return alpha_6 parameter from Bandt and Shiha's Order Patterns in Time 
    Series paper for fractional Brownian motion.
    
    Keyword arguments:
    H -- Hurst parameter H; H = 1/2 is ordinary brownian motion

        >>> _alpha6(0.5)
        -0.5773502691896258

    """

    return (2**(2*H) - 3**(2*H) - 1) / (2*3**H)

def _alpha7(H):
    """
    Return alpha_7 parameter from Bandt and Shiha's Order Patterns in Time 
    Series paper for fractional Brownian motion.
    
    Keyword arguments:
    H -- Hurst parameter H; H = 1/2 is ordinary brownian motion

        >>> _alpha7(0.5)
        0.0

    """

    return (3**(2*H) - 2**(2*H) - 1) / 2**(H+1)

def _alpha8(H):
    """
    Return alpha_8 parameter from Bandt and Shiha's Order Patterns in Time 
    Series paper for fractional Brownian motion.
    
    Keyword arguments:
    H -- Hurst parameter H; H = 1/2 is ordinary brownian motion

        >>> _alpha8(0.5)
        0.5773502691896258

    """

    return (2**(2*H) - 1) / 3**H

def _reverse(permutation):
    """
    Return the reverse of a permutation given in 1-line notation.

    Keyword arguments:
    permutation -- A permutation's 1-line notation

        >>> _reverse((2,4,3,1))
        (1, 3, 4, 2)

    """

    return permutation[::-1]

def _complement(permutation):
    """
    Return the complement of a permutation of {1,...,n}

    Keywork arguments:
    permutation -- A permutation of {1,...,n} in 1-line notation.

        >>> _complement((2,4,3,1))
        (3, 1, 2, 4)

    """
    return tuple(len(permutation) + 1 - x for x in permutation)

def _reverse_complement(permutation):
    """
    Return the reverse complement of a permutation of {1,...,n}

    Keyword arguments:
    permutation -- A permutation of {1,...,n}

        >>> _reverse_complement((2,4,3,1))
        (4, 2, 1, 3)

    """

    return _reverse(_complement(permutation))

def _gaussian4_dictionary(H):
    """
    Return a dictionary comprised of ordinal pattern probabilities for a
    fractional Brownian motion based on Hurst parameter H. The dictionary
    does not contain all permutations of S_4. One must apply reverse,
    complement, and reverse complement to get the rest.

    Keyword arguments:
    H -- Hurst parameter for the fractional Brownian motion

        >>> _gaussian4_dictionary(0.5)
        {(1, 2, 3, 4): 0.125, (3, 1, 4, 2): 0.014623304674318438, (4, 2, 3, 1): 0.04166666666666666, (2, 1, 4, 3): 0.02704336199234815, (1, 2, 4, 3): 0.062499999999999986, (1, 4, 2, 3): 0.0208333333333333, (3, 1, 2, 4): 0.035456638007651795, (1, 4, 3, 2): 0.02704336199234815}

    """

    p = {}

    p[(1,2,3,4)] = 1/8 + 1/(4*pi) *\
                   (asin(_alpha1(H)) + 2*asin(_alpha2(H)))
    p[(3,1,4,2)] = 1/8 + 1/(4*pi) *\
                   (2*asin(_alpha3(H)) + asin(_alpha4(H)))
    p[(4,2,3,1)] = 1/8 + 1/(4*pi) *\
                   (asin(_alpha4(H)) - 2*asin(_alpha5(H)))
    p[(2,1,4,3)] = 1/8 + 1/(4*pi) *\
                   (2*asin(_alpha6(H)) + asin(_alpha1(H)))
    p[(1,2,4,3)] = 1/8 + 1/(4*pi) *\
                   (asin(_alpha7(H)) - asin(_alpha1(H)) - asin(_alpha5(H)))
    p[(1,4,2,3)] = 1/8 + 1/(4*pi) *\
                   (asin(_alpha7(H)) - asin(_alpha4(H)) - asin(_alpha5(H)))
    p[(3,1,2,4)] = 1/8 + 1/(4*pi) *\
                   (asin(_alpha3(H)) + asin(_alpha8(H)) - asin(_alpha5(H)))
    p[(1,4,3,2)] = 1/8 + 1/(4*pi) *\
                   (asin(_alpha6(H)) - asin(_alpha8(H)) + asin(_alpha2(H)))

    return p
