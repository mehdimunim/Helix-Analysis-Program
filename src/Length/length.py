from numpy.linalg import norm


def length(axis, alpha_carbons):
    """
    Calculate the length of the structure

    Can serve to learn about compression and extension of a helix during MD simulation
    cf. TRAJELIX 
    """
    _, first_orig = normal(axis, alpha_carbons[0])
    _, last_orig = normal(axis, alpha_carbons[-1])

    return distance(first_orig, last_orig)
