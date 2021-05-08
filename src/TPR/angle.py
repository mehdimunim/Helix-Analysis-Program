from numpy import dot
from numpy.linalg import norm
from math import acos


def angle(v1, v2):
    """
    Calculate the angle between v1 and v2.
    ---
    Return:
    angle: angle in rad in [0, pi]

    """
    cos_angle = dot(v1, v2)/(norm(v1)*norm(v2))
    angle = acos(cos_angle)
    return angle


def test_angle():
    from math import pi
    import numpy as np

    t = 10**(-5)

    v1 = np.array([1, 2, 3])
    v2 = np.array([-1, -2, -3])

    ang1 = angle(v1, v2)
    print("pi expected: ", ang1)

    ang2 = angle(v1, v1)
    print("0 expected: ", ang2)


test_angle()
