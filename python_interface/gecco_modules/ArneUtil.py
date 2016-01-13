def combine_dicts(first,second):
    """concatenates two dictionaries

    if both have an identical key. the first dictionary takes precedence
    """

    def _combine_dicts_core(alpha,beta,x):
        if x in alpha:
            return alpha.get(x)
        else:
            return beta.get(x)
    return {x: _combine_dicts_core(first,second,x)  for x in set(first).union(second)}




def ggT(a,b):
    """ calculates the largest common denominator"""
    if b == 0 :
        return a
    else:
        return (b,a%b)
