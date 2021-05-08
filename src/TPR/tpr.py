#
# Turn per residue
#

from normal import normal
from angle import angle
from fit import fit


def tpr(alpha_carbons, axis_direction, axis_center):
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
    axis_direction: inertia axis direction of the structure (often an alpha-helix)
    axis_center: center of the intertia axis

    ---
    Return:
    theta : turn per angle per residue (in rad)

    """

    # First step
    # Calculate the turn angles between first residue and i

    phi = []
    phi_i = 0
    for i in range(1, len(alpha_carbons)):
        # calculate normal vectors for residue i-1 and i
        vec_before, _ = normal(axis_direction, axis_center, alpha_carbons[i-1])
        vec_after, _ = normal(axis_direction, axis_center, alpha_carbons[i])

        # angle between i-1 and i
        angle_between = angle(vec_before, vec_after)

        # angle between first residue and i
        phi_i += angle_between

        phi.append(phi_i)

    # Second step
    # Find turn angle per residue with linear regression

    theta = fit(phi)

    return theta
