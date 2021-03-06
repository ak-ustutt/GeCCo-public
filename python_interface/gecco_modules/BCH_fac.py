import math

def set_BCH_factor(n, k, use_minus = True):
    """Return the factor for the BHC expansion
    """
    fac = 1.0/(math.factorial(k)*math.factorial(n-k))
    if (use_minus and k % 2 == 1):
        fac = -fac
    return fac
