import numpy as np


def fit(phi):
    """
    Find best phi_0 and theta such as:

    phi_i = phi_0 + i*theta

    ---
    Parameters:
    phi: list of angles between first residue and i

    ---
    Return:
    theta: the turn angle per residue 

    """
    x = [i for i in range(len(phi))]

    _, theta = np.polynomial.polynomial.Polynomial.fit(
        x, phi, deg=1).convert().coef

    return theta
