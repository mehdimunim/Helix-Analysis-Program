from dssp import *
from axis import *
from tpr import *
from parse import parse
import math


def test_tpr():
    backbone = parse("glut1.pdb")
    dssp = DSSP(backbone)
    list_helices = dssp.get_ca()

    for helix in list_helices:
        orig, axis = principal_axis(helix)
        theta = tpr(helix, axis, orig)
        # The angle is expected to be about 100° as there are 3.6 residue by alpha-helix turn
        print("theta = {:.1f}°".format(theta*180/math.pi))


test_tpr()
