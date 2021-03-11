import math


def radcirc(a, b, c):
    """

    ***
    Parameters
    a, b, c: three vectors

    Returns:
    r: the radius of the circle going through a,b and c
    """
    rHB, rb, rac, angabc = angdistw(b, a, c)
    sabc = math.sin(angabc)
    if (sabc < 0.00001):
        r = 999999.9
    else:
        r = rac/(math.sin(angabc)*2)
    return r
