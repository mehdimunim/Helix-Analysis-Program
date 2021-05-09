from dssp import *
from axis import *
from tpr import *
import math


def test_tpr():
    dssp = DSSP("glut1.pdb")
    list_helices = dssp.get_ca()

    for helix in list_helices:
        orig, axis = principal_axis(helix)
        theta = tpr(helix, axis, orig)
        # The angle is expected to be about 100° as there are 3.6 residue by alpha-helix turn
        print("theta = {:.1f}°".format(theta*180/math.pi))


test_tpr()
