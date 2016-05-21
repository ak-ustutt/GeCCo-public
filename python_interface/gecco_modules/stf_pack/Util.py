def combine_dicts(first,second):
    """concatenates two dictionaries

    if both have an identical key. the first dictionary takes precedence
    """
    combined = {}
    combined.update(second)
    combined.update(first)
    return combined




def ggT(a,b):
    """ calculates the largest common denominator

    uses Euclids (?) Algorithm
    """
    if b == 0 :
        return a
    else:
        return ggT(b,a%b)
