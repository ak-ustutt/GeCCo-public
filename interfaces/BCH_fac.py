import math
# The factor of BCH expansion
def set_BCH_factor(n, k, use_minus = True):
    fac = 1.0/(math.factorial(k)*math.factorial(n-k))
    if (use_minus and k % 2 == 1):
        fac = -fac
    return fac
