#
# Turn per residue
#

from normal import normal
from angle import angle
from regression import regression


def tpr(alpha_carbons, axis):
    """
    Calculate the turn angle per residue;
    Inspired from TRAJELIX from simulaid.

    Two steps:

    1- Calculate the turn angles phi(i) between first residue and i

    2- Fit ph(i) = theta*i + phi_0 where theta is the turn angle per residue

    (Mezei, Filizola, 2006)

    ---
    Parameters:
    alpha_carbons: alpha carbons (used to calculate the angle)
    axis: inertia axis of the structure (often an alpha-helix)

    ---
    Return:
    theta : turn per angle per residue

    """

    # First step
    # Calculate the turn angles between first residue and i

    phi = []
    phi_i = 0
    for i in range(3, n):
        # calculate normal vectors for residue i-1 and i
        vec_before = normal(axis, alpha_carbons[i-1])
        vec_after = normal(axis, alpha_carbons[i])

        # angle between i-1 and i
        angle = angle(vec_before, vec_after)

        # angle between first residue and i
        phi_i += angle

        phi.append(phi_i)

    # Second step
    # Find turn angle per residue with linear regression

    theta, phi_0 = linear_regression(phi)

    return theta
