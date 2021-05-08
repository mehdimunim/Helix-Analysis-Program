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
        print("theta = ", theta*180/math.pi)


test_tpr()
