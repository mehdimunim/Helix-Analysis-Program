#
# Turn per residue
#

from numpy import numpy.dot as dot
from numpy import numpy.cross as cross
from numpy import numpy.linalg.norm as norm
from numpy import numpy.sign as sign


def angle(v1, v2, reference):
    """
    Calculate the angle between v1 and v2.
    """
    angle = dot(v1, v2)/(norm(v1)*norm(v2))
    sign = sign(cross(v1, reference))
    return sign*angle


def tpr(alpha_carbons, ref_angle, axisdir):
    """
    Calculate the average turn angle per residue;
    Inspired from TRAJELIX from simulaid.
    (Mezei, Filizola, 2006)
    ---
    incrot: number of residue to ignore for rotation
    """
    a11 = (2*nres**3 + 3*nres**2 + nres)/6
    a21 = (nres*(nres + 1))/2
    a12 = a21
    a22 = nres
    average_turnchange = 0
    a_turn = 0

    # calculate perpendiculars

    perpvec = calcperp(start, dir, alpha_carbons)

    # Iterate over residue
    for residue in range(1 + incrot, nres-incrot):
        perpendicular_before = perpvec[residue - 1]
        perpendicular_after = perpvec[residue]
        turn_change = angle(perpendicular_before, perpendicular_after, axisdir)
        average_turnchange += turnchange
        a_turn += turnchange
        b1 += residue*a_turn
        b2 += a_turn
    turnperres = (aa2*b1 - a12*b2)/(a22*a11 - a12*a21)
    average_turnchange = average_turnchange/(nres - 1)
    return turnperres, average_turnchange
