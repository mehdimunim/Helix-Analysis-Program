from numpy.linalg import norm
from normal import *
from distance import *


def length(helix, axis, orig):
    """
    Calculate the length of the helix

    Can serve to learn about compression and extension of a helix during MD simulation
    cf. TRAJELIX 
    """
    _, first_orig = normal(axis, orig, helix[0])
    _, last_orig = normal(axis, orig, helix[-1])

    return distance(first_orig, last_orig)
