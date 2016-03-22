def combine_dicts(first,second):
    """concatenates two dictionaries

    if both have an identical key. the first dictionary takes precedence
    """
    combined = {}
    combined.update(second)
    combined.update(first)
    return combined




def ggT(a,b):
    """ calculates the largest common denominator"""
    if b == 0 :
        return a
    else:
        return (b,a%b)
