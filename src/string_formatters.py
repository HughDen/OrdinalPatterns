def float_dictionary_string(D):
    r"""
    Return a nicer-than-default string from a dictionary that maps to floats.

    Keyword arguments:
    D -- A dictionary whose keys can be converted to strings and whose values
         are floats

        >>> float_dictionary_string({(1,2): 0.5, (2,1): 0.5})       
        '(1, 2) -> 0.5000     (2, 1) -> 0.5000     \n'
    """
    counter = 0
    s = ""
    for k in D.keys():
        if counter % 3 == 2:
            s += "%s -> %.4f \n" % (str(k), D[k])
        else:
            s += "%s -> %.4f     " % (str(k), D[k])
        counter += 1
        
    if len(D) < 3:
        s += "\n"

    return s

def int_dictionary_string(D):
    r"""
    Return a nicer-than-default string from a dictionary that maps to ints.

    Keyword arguments:
    D -- A dictionary whose keys can be converted to strings and whose values
         are ints

        >>> int_dictionary_string({(1,2): 1, (2,1): 1})       
        '(1, 2) -> 1     (2, 1) -> 1     \n'
    """
    counter = 0
    s = ""
    for k in D.keys():
        if counter % 3 == 2:
            s += "%s -> %d \n" % (str(k), D[k])
        else:
            s += "%s -> %d     " % (str(k), D[k])
        counter += 1
        
    if len(D) < 3:
        s += "\n"

    return s

def matrix_list_string(matrix_list):
    r"""
    Return a nicer-than-default string from a list of matrices.

    Keyword arguments:
    matrix_list -- a list of matrices with the same number of rows

        >>> matrix_list_string([[[7,1],[3,2]],[[5,4],[6,-1]]])       
        '7  1    5  4    \n3  2    6  -1    \n'
    """
    SPACING = 3
    s = ""
    N = len(matrix_list)
    rounds = N // SPACING

    for q in range(rounds):
        # This might fail badly if the matrices have different row lengths
        for row_number in range(len(matrix_list[SPACING*q])):
            for r in range(SPACING):
                s += '  '.join([str(x) for x in matrix_list[SPACING*q + r][row_number]])
                s += '    '
            s += "\n"
        s += "\n"

    if N > rounds*SPACING:
        for row_number in range(len(matrix_list[SPACING*rounds])):
            for r in range(rounds*SPACING, N):
                s += '  '.join([str(x) for x in matrix_list[r][row_number]])
                s += '    '
            s += '\n'
           
    return s

def address_list_string(address_list):
    r"""
    Return a string for a list of addresses as a string of representing 
    matrices.

    Keyword arguments:
    address_list -- a list of addresses of alcoves

        >>> address_list_string([{(1,2):0, (1,3):0, (2,3):0}, {(1,2):0, (1,3):1, (2,3):0}])       
        '0  0  0    0  0  1    \n0  0  0    0  0  0    \n0  0  0    0  0  0    \n'
    
    """

    matrix_list = [_address_to_matrix(address) for address in address_list]
    return matrix_list_string(matrix_list)

def _address_to_matrix(address):
    """
    Return a matrix associated to the address for an alcove. Addresses
    are maps from pairs to values. Thus, they may be represented in matrix
    form.

    Keyword arguments:
    address -- the address for an alcove.

        >>> _address_to_matrix({(1,2):0,(1,3):1,(2,3):0})
        [[0, 0, 1], [0, 0, 0], [0, 0, 0]]

    """

    K = [x for k in address.keys() for x in k]
    N = max(K)

    A = [[0 for _ in range(N)] for _ in range(N)]

    for k in address.keys():
        i,j = k
        A[i-1][j-1] = address[k]

    return A
