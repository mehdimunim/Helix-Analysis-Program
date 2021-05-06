#
# Turn per residue
#

from normal import normal
from angle import angle
from fit_angles import fit_angles


def tpr(alpha_carbons, axis):
    """
    Calculate the turn angle per residue;
    Inspired from TRAJELIX from simulaid.

    Two steps:

    1- Calculate the turn angles phi_i between first residue and i

    2- Fit phi_i = theta*i + phi_0 where theta is the turn angle per residue

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
    for pos in range(1, len(alpha_carbons)):
        # calculate normal vectors for residue i-1 and i
        vec_before, _ = normal(axis, alpha_carbons[i-1])
        vec_after, _ = normal(axis, alpha_carbons[i])

        # angle between i-1 and i
        angle = angle(vec_before, vec_after)

        # angle between first residue and i
        phi_i += angle

        phi.append(phi_i)

    # Second step
    # Find turn angle per residue with linear regression

    theta, phi_0 = fit_angles(phi)

    return theta
