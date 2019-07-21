from itertools import permutations, combinations, product
from collections import defaultdict, deque
from math import floor, ceil, sqrt, pi, gamma, factorial

def hypercube_address_generator(n):
    '''
    Return a generator for addresses of type A alcoves in [0,1]^n.

    Keyword arguments:
    n -- Dimension of [0,1]^n

        >>> tuple(x for x in hypercube_address_generator(2))
        (((0, 0), (1,)), ((0, 0), (0,)))

    '''

    for p in permutations(range(1,n+1)):
        yield address_from_permutation(p)

def address_dictionary_generator(n):
    '''
    Return a generator for addresses of type A alcoves in [0,1]^n
    represented as dictionaries with roots (i,j) as keys.

    Keyword arguments:
    n -- Dimension of [0,1]^n

        >>> tuple(x for x in address_dictionary_generator(2))
        ({(1, 2): 0, (1, 3): 1, (2, 3): 0}, {(1, 2): 0, (1, 3): 0, (2, 3): 0})

    '''

    for p in permutations(range(1,n+1)):
        yield dictionary_of_address(address_from_permutation(p))

def BFS_weak_order_addresses(permutation):
    '''
    Return a tuple of addresses (represented as dictionaries) for the 
    polytope P_pi constructed from a given permutation.

    Keyword arguments:
    permutation -- a permutation of 1,2,...,n

        >>> BFS_weak_order_addresses((1,3,2))
        ({(1, 2): 0, (1, 3): 0, (2, 3): 0},)

    '''

    N = len(permutation)
    weak_order_dictionaries = []
    ideal = root_ideal(permutation)
    base_dictionary = {x : 0 for x in ideal}

    Q = deque([(base_dictionary,(1,3))])

    while Q:
        dictionary, current_root = Q.popleft()

        # If all roots have values, it is one of the addresses in the interval
        if len(dictionary.keys()) == (N*(N-1)) // 2:
            weak_order_dictionaries.append(dictionary)
        else:
            x,y = current_root
            
            # All possibilities require the next root; let's do that first
            # We move along all pairs with the same difference (i.e. height)
            # until we reach the end. At the end, the difference goes up and
            # we start back at (1,_)
            
            if y == len(permutation):
                new_delta = y - x + 1
                new_x = 1
                new_y = 1 + new_delta
            else:
                new_x = x + 1
                new_y = y + 1
            
            # If the current root is not there, we add all
            # possibilities satisfying the constraints of Shi's Theorem
            
            if current_root not in dictionary.keys():
                sums = set([dictionary[(x,a)] + dictionary[(a,y)]
                        for a in range(x+1,y)])
                
                
                
                if len(sums) == 2:
                    dictionary[(x,y)] = max(sums)
                    Q.append((dictionary,(new_x,new_y)))
                elif len(sums) == 1:
                    new_dictionary = {k:dictionary[k]
                                      for k in dictionary.keys()}
                    new_dictionary[(x,y)] = max(sums) + 1
                    dictionary[(x,y)] = max(sums)
                    Q.append((dictionary,(new_x,new_y)))
                    Q.append((new_dictionary,(new_x,new_y)))
                else:
                    print("This is not happenning")
            else:
                Q.append((dictionary,(new_x,new_y)))
                
    return tuple(weak_order_dictionaries)

def root_ideal(permutation):
    '''
    Return the ideal of consecutive roots of a permutation represented as pairs.

    Keyword arguments:
    permutation -- a permutation of 1,2,...,n
    
        >>> root_ideal((1,3,2,4))
        ((1, 2), (1, 3), (2, 3), (3, 4), (2, 4))
    '''
    
    I = set([])
    
    for i in range(len(permutation) - 1):
        low = min(permutation[i],permutation[i+1])
        high = max(permutation[i],permutation[i+1])
        # Add all pairs in range determined by consecutive values
        for j in range(low, high):
            for k in range(j+1,high+1):
                I.add((j,k))

    return tuple(I)
    
# A consequence of Shi's Theorem and the translational symmetry of the affine
# hyperplane arrangement is that any address outside of [0,1]^n
# can be obtained from an address in [0,1]^n by a suitable translation. The
# method region_alcove_generator utilizes this to generate alcoves in [a,b]^n.

def region_alcove_generator(n, low, high):
    '''
    Return a generator for all alcoves whose address on simple roots is 
    in [low,high]^n.

    Keyword arguments:
    n -- positive integer for dimension of [0,1]^n
    low -- integer for left endpoint of [low,high]
    high -- integer for right endpoint of [low,high]

        >>> tuple(region_alcove_generator(2, -1, 0))
        (((-1, -1), (-1,)), ((-1, -1), (-2,)), ((-1, 0), (0,)), ((-1, 0), (-1,)), ((0, -1), (0,)), ((0, -1), (-1,)), ((0, 0), (1,)), ((0, 0), (0,)))

    '''
    
    for base in product(range(low,high + 1), repeat = n):
        base_alcove = _translation_base(base)
        
        for delta_alcove in hypercube_address_generator(n):
            yield alcove_add(delta_alcove, base_alcove)


def alcove_add(alcove1, alcove2):
    '''
    Return entry-wise addition of the addresses of two alcoves in tuple form. 
    
    Keyword arguments:
    alcove1 -- The address of an alcove represented as a tuple of tuples
    alcove2 -- The address of an alcove represented as a tuple of tuples

        >>> alcove_add(((-1,-1),(-2,)), ((0,0),(1,)))
        ((-1, -1), (-1,))

    '''
    return tuple(tuple(alcove1[i][j] + alcove2[i][j]
                         for j in range(len(alcove1) - i))
                  for i in range(len(alcove1)))

# The method address_from_permutation is based on Richard Stanley's 'Eulerian 
# partitions of a hypercube'. In that paper, he gives a volume-preserving
# transformation that places a sum y_1 + ... + y_k between k and k + 1, where
# k is the number of ascents (including the 'ascent by convention' at 0).
#
# Our main modification is to include hyperplanes of the form y_i + ... + y_j.
# This gives a bijection between addresses in [0,1]^n and permutations. It's 
# straightforward to implement and reasonably fast.
#
# To get a dictionary from positive roots as pairs (i,j) to nonnegative
# integers, use the method dictionary_of_addresses.

def address_from_permutation(permutation):
    '''
    Return the address of an alcove in [0,1]^n corresponding to the given 
    permutation. 

    Keyword arguments
    permutation -- A tuple or list of comparable items

        >>> address_from_permutation((1,3,2))
        ((0, 0, 0), (1, 0), (1,))
    '''
    n = len(permutation)
    ac = _cumulative_ascent_count(permutation)
    mod_perm = (0,) + permutation
    
    return tuple(tuple((ac[i+j] - ac[j]) - (mod_perm[j] < mod_perm[i+j])
                       for j in range(n+1 - i)) for i in range(1,n+1))

def _cumulative_ascent_count(permutation):
    '''
    Return a list of the number of descents of a permutation up to position i.

    Keyword arguments:
    permutation -- A tuple or list of comparable items

        >>> _cumulative_ascent_count((3,5,1,2,4))
        (0, 1, 2, 2, 3, 4)
    '''

    ascent_count = [0]
    mod_perm = (0,) + permutation
    
    count = 0
    for i in range(len(mod_perm) - 1):
        if mod_perm[i] < mod_perm[i+1]:
            count += 1
        ascent_count.append(count)

    return tuple(ascent_count)

def dictionary_of_address(address):
    '''
    Return a dictionary representing the given address as a dictionary whose
    keys are pairs (i,j) that represent positive roots. 

    Keyword arguments:
    address -- A tuple of tuples representing an alcove's address

        >>> dictionary_of_address(((0,0,0),(1,0),(1,)))
        {(1, 2): 0, (1, 3): 1, (1, 4): 1, (2, 3): 0, (2, 4): 0, (3, 4): 0}

    '''

    address_dictionary = {}
    for (i,j) in combinations(range(1, len(address)+2), 2):
        a = j - i - 1
        b = i - 1
        address_dictionary[(i,j)] = address[a][b]

    return address_dictionary

def address_to_ordinal_pattern(address):
    '''
    Return the ordinal pattern associated to the given address.

    Keyword arguments:
    address -- A tuple of tuples representing an alcove's address

        >>> address_to_ordinal_pattern(((-1,0,0),(-1,0),(0,)))
        (3, 1, 2, 4)

    '''
    
    n = len(address)
    code = [0 for _ in range(n+1)]

    for i in range(n):
        for j in range(n - i):
            if address[i][j] < 0:
                code[j] += 1

    return code_to_permutation(code)

def code_to_permutation(code):
    '''
    Return the permutation for the given code. The code here is the standard
    Lehmer code in which each entry records the number of inversions from
    that position.

    Keyword arguments:
    code -- A tuple or list from [0,n-1] X [0,n-2] X ... {0}

        >>> code_to_permutation((2,0,0,0))
        (3, 1, 2, 4)

    '''

    values = list(range(1,len(code) + 1))
    permutation = []

    for i in range(len(code)):
        permutation.append(values[code[i]])
        del values[code[i]]

    return tuple(permutation)

def gaussian_estimate(n, R):
    """
    Return a dictionary mapping permutations of S_n to an 
    estimate of the ordinal pattern probability of a Gaussian
    walk.

    Keyword arguments:
    n -- dimension of the n-sphere
    R -- radius of the n-sphere

        >>> gaussian_estimate(2,8)
        {(2, 1): 0.5000000000000001, (1, 2): 0.5000000000000001}

    """

    count = _gaussian_count(n, R)
    volume = _hypersphere_volume(n-1,R)

    highs = {p: count[p][1] / (volume*factorial(n-1))
             for p in count.keys()}

    lows = {p:count[p][0] / (volume*factorial(n-1))
            for p in count.keys()}

    highsums = sum(highs.values())
    lowsums = sum(lows.values())

    if highsums != lowsums:
        weight = (highsums - 1) / (highsums - lowsums)
    else:
        weight = 0.5
        
    estimates = {p: (weight*count[p][0] + (1-weight)* count[p][1]) / (volume*factorial(n-1))
                 for p in count.keys()}

    return {p: estimates[p] for p in count.keys()}
    
def _gaussian_count(n, R):
    """
    Return a dictionary mapping permutations of S_n to counts
    that underestimate and overestimate volumes of D_pi 
    inside an n-sphere of radius R.

    Keyword arguments:
    n -- dimension of n-sphere
    R -- radius of n-sphere

        >>> _gaussian_count(3, 8)
        {(3, 2, 1): (82, 112), (3, 1, 2): (41, 56), (2, 1, 3): (41, 56), (2, 3, 1): (41, 56), (1, 3, 2): (41, 56), (1, 2, 3): (82, 112)}

    """

    undercount = defaultdict(int)
    overcount = defaultdict(int)
    
    for base in integer_sphere_generator(n-1, R):
        base_alcove = _translation_base(base)
        
        for delta_alcove in hypercube_address_generator(n-1):
            alcove = alcove_add(base_alcove, delta_alcove)
            overcount[address_to_ordinal_pattern(alcove)] += 1
            if sum([(b + (b >= 0)) ** 2 for b in base]) <= R**2:
                undercount[address_to_ordinal_pattern(alcove)] += 1

    return {p:(undercount[p],overcount[p]) for p in overcount.keys()}
            
def integer_sphere_generator(n,R):
    """
    Return a generator of integer n-tuples such that the sum of 
    squares of the minimum of values x and x + 1 in the tuple are
    less than R**2.

    Keyword arguments:
    n -- dimension of the generated tuples
    R -- the radius of an n-sphere

        >>> tuple(integer_sphere_generator(2, 1.35))
        ((-2, -1), (-2, 0), (-1, -2), (-1, -1), (-1, 0), (-1, 1), (0, -2), (0, -1), (0, 0), (0, 1), (1, -1), (1, 0))

    """
    
    if n == 1:
        for x in range(floor(-R), ceil(R)):
            yield (x,)
    else:
        for x in range(floor(-R), ceil(R)):
            if x >= 0:
                for v in integer_sphere_generator(n - 1, sqrt(R**2 - x**2)):
                    yield (x,) + v
            else:
                for v in integer_sphere_generator(n - 1, sqrt(R**2 - (x + 1)**2)):
                    yield (x,) + v


def _hypersphere_volume(n,R):
    """
    Return the volume of an n-dimensional hypersphere of radius R.

    Keyword arguments:
    n -- dimension of the hypersphere
    R -- radius of the hypersphere

        >>> _hypersphere_volume(3,4)
        268.082573106329

    """

    return R**n * pi**(n/2) / gamma(1 + n / 2)

# For fixed values of address values on simple roots, the partial sums
# determine the minimum values of all positive roots.

def _translation_base(address_values):
    """
    Return an address (represented by a tuple of tuples) that has the
    smallest values on roots for a given tuple of address values on
    simple roots.

    Keyword arguments
    address_values -- a tuple containing values of an address on 
                      simple roots.

        >>> _translation_base((2,-3,4,-1))
        ((2, -3, 4, -1), (-1, 1, 3), (3, 0), (2,))

    """

    sums = _simple_partial_sums(address_values)
    return tuple(tuple(sums[i+j] - sums[j]
                       for j in range(len(sums)-i))
                 for i in range(1,len(sums)))
    
def _simple_partial_sums(address_values):
    """
    Return a tuple whose i-th entry is the i-th partial sum of 
    address values on simple roots.

    Keyword arguments:
    address_values -- a tuple containing values of an address on
                      simple roots
        
        >>> _simple_partial_sums((2,-3,4,-1))
        (0, 2, -1, 3, 2)

    """

    return tuple(sum(address_values[:i])
                 for i in range(len(address_values) + 1))
