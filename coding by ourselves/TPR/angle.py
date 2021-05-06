from numpy import dot
from numpy.linalg import norm
from math import arccos


def angle(v1, v2):
    """
    Calculate the angle between v1 and v2.
    """
    cos_angle = dot(v1, v2)/(norm(v1)*norm(v2))
    angle = arccos(cos_angle)
    return angle
