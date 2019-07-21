import random
import numpy as np
from functools import partial
from collections import defaultdict
from itertools import permutations

def ordinal_laplace_walk(n, iterations):
    """
    Return a dictionary that maps permutations in S_n to their
    frequency of occurrence in a random walk simulation with a specified
    number of iterations where the steps are drawn from a symmetric uniform
    distribution.

    Keyword arguments:
    n -- length of the ordinal patterns
    iterations -- number of simulations to run

        >>> import random
        >>> random.seed(27) 
        >>> ordinal_laplace_walk(3, 25000)
        {(1, 2, 3): 0.25428, (1, 3, 2): 0.12364, (2, 1, 3): 0.122, (2, 3, 1): 0.12684, (3, 1, 2): 0.12572, (3, 2, 1): 0.24752}

    """

    # A random variable drawn from the Laplace distribution is 
    # the same as the subtraction of two exponential random variables.

    expo1 = partial(random.expovariate, 1)
    expo2 = partial(random.expovariate, 1)

    laplace = lambda : expo1() - expo2()
    return ordinal_simulation_walk(n, laplace, iterations)

def ordinal_uniform_walk(n, iterations):
    """
    Return a dictionary that maps permutations in S_n to their
    frequency of occurrence in a random walk simulation with a specified
    number of iterations where the steps are drawn from a symmetric uniform
    distribution.

    Keyword arguments:
    n -- length of the ordinal patterns
    iterations -- number of simulations to run

        >>> import random
        >>> random.seed(8) 
        >>> ordinal_uniform_walk(3, 25000)
        {(1, 2, 3): 0.2452, (1, 3, 2): 0.12284, (2, 1, 3): 0.12944, (2, 3, 1): 0.12404, (3, 1, 2): 0.12272, (3, 2, 1): 0.25576}

    """
    uniform = partial(random.uniform, -1, 1)
    return ordinal_simulation_walk(n, uniform, iterations)

def ordinal_gaussian_walk(n, iterations):
    """
    Return a dictionary that maps permutations in S_n to their
    frequency of occurrence in a random walk simulation with a specified
    number of iterations where the steps are drawn from a symmetric Gaussian.

    Keyword arguments:
    n -- length of the ordinal patterns
    iterations -- number of simulations to run

        >>> import random
        >>> random.seed(17)
        >>> ordinal_gaussian_walk(3, 25000)
        {(1, 2, 3): 0.25412, (1, 3, 2): 0.12296, (2, 1, 3): 0.12164, (2, 3, 1): 0.12564, (3, 1, 2): 0.12348, (3, 2, 1): 0.25216}

    """

    gauss = partial(random.gauss, 0, 1)
    return ordinal_simulation_walk(n, gauss, iterations)

def ordinal_simulation_walk(n, generator, iterations):
    """
    Return a dictionary that maps permutations in S_n to their
    frequency of occurrence in a random walk simulation with a specified 
    number of iterations from a user-supplied pseudo-random number generator.

    Keyword arguments:
    n -- length of the ordinal patterns
    generator -- pseduo-random number generator
    iterations -- number of simulations to run

    For standard tools like random.gauss, the generality of passing a random
    number generator means passing a function with no parameters, as
    shown below.
 
        >>> import random
        >>> from functools import partial
        >>> gauss = partial(random.gauss, 0, 1)
        >>> uniform = partial(random.uniform, -1, 1)
        >>> random.seed(17)
        >>> ordinal_simulation_walk(3, gauss, 25000)
        {(1, 2, 3): 0.25412, (1, 3, 2): 0.12296, (2, 1, 3): 0.12164, (2, 3, 1): 0.12564, (3, 1, 2): 0.12348, (3, 2, 1): 0.25216}
        >>> random.seed(8) 
        >>> ordinal_simulation_walk(3, uniform, 25000)
        {(1, 2, 3): 0.2452, (1, 3, 2): 0.12284, (2, 1, 3): 0.12944, (2, 3, 1): 0.12404, (3, 1, 2): 0.12272, (3, 2, 1): 0.25576}
 
    """

    pattern_count = defaultdict(int) #defaultdict ensures default values are 0.

    for _ in range(iterations):
        deltalist = [generator() for _ in range(n - 1)]
        floatlist = [sum(deltalist[:i]) for i in range(n)]
        pattern_count[ordinal_pattern(floatlist)] += 1

    return {p:(pattern_count[p] / iterations) 
            for p in permutations(range(1,n+1))}

def ordinal_simulation(n, generator, iterations):
    """
    Return a dictionary that maps permutations in S_n to their
    frequency of occurrence in a simulation with a specified number of
    iterations from a user supplied pseudo-random number generator.

    Keyword arguments:
    n -- length of the ordinal patterns
    generator -- pseduo-random number generator
    iterations -- number of simulations to run

    For standard tools like random.gauss, the generality of passing a random
    number generator means passing a function with no parameters, as
    shown below.
 
        >>> import random
        >>> from functools import partial
        >>> gauss = partial(random.gauss, 0, 1)
        >>> exponential = partial(random.expovariate, 1)
        >>> random.seed(17)
        >>> ordinal_simulation(2, gauss, 25000)
        {(1, 2): 0.4988, (2, 1): 0.5012}
        >>> random.seed(8)
        >>> ordinal_simulation(3, exponential, 25000)
        {(1, 2, 3): 0.16856, (1, 3, 2): 0.16792, (2, 1, 3): 0.16816, (2, 3, 1): 0.16616, (3, 1, 2): 0.16332, (3, 2, 1): 0.16588}
 
    """

    pattern_count = defaultdict(int) #defaultdict ensures default values are 0.

    for _ in range(iterations):
        floatlist = [generator() for _ in range(n)]
        pattern_count[ordinal_pattern(floatlist)] += 1

    return {p:(pattern_count[p] / iterations) 
            for p in permutations(range(1,n+1))}

def ordinal_pattern(floatlist):
    """ 
    Return the ordinal pattern from a list of comparable elements.
    Our intended use-case is a list of distinct floats. The routine returns
    a sensible answer for other comparable data types as well.

    The 1-line representation of the permutation is returned as a tuple of 
    ints in {1,...,n}.

    Keyword argument:
    floatlist -- a list of floats to convert to its ordinal pattern   

        >>> ordinal_pattern([-1.2, 2.7, 1.5])
        (1, 3, 2)
        >>> ordinal_pattern([8000,4000,1000,2000])
        (4, 3, 1, 2)
       
    """

    # We get indices that could be used to sort the original array from
    # the argsort() method. The inverse permutation of those indices 
    # is obtained by argsort(), which produces the original order.

    sorting_indices = np.argsort(floatlist)
    inverse_of_indices = np.argsort(sorting_indices)

    return tuple(x + 1 for x in inverse_of_indices)    
