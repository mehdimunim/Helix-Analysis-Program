from length import *
from axis import *
from dssp import *


def test_length():
    dssp = DSSP("glut1.pdb")
    list_helices = dssp.get_ca()
    for helix in list_helices:
        orig, axis = principal_axis(helix)
        len = length(helix, axis, orig)
        print(len)


test_length()
