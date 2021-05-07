import math


def distance(atom1, atom2):
    """
    Calculates distance between atom1 and atom2
    atom1, atom2 are 3-arrays of coordinates
    """
    return math.sqrt((atom1[0] - atom2[0])**2 + (atom1[1] - atom2[1])**2 + (atom1[2] - atom2[2])**2)
