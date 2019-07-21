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
        
    if len(D) < 6:
        s += "\n"

    return s
