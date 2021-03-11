import numpy as np


def normplane(a, b, c):
    """
    ***
    Parameters:
    a,b,c: three vectors

    Returns:
    rn the normal to the plane of a,b and c
    """
    d1 = a - b
    rn = np.cross(d1, d2)
    rn = rn/np.linalg.norm(rn)
    return rn
