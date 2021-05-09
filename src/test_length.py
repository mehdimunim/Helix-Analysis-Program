from length import *
from axis import *
from dssp import *
from parse import parse


def test_length():
    import os
    os.chdir(
        "C:\\Users\\Mehdi\\Documents\\GitHub\\Interdisciplinary-Project\\resource")
    backbone = parse("glut1.pdb")
    dssp = DSSP(backbone)
    list_helices = dssp.get_ca()
    for helix in list_helices:
        orig, axis = principal_axis(helix)
        len = length(helix, axis, orig)
        print(len)

def test_length_traj():
    import os
    os.chdir(
        "C:\\Users\\Mehdi\\Documents\\GitHub\\Interdisciplinary-Project\\resource")
    backbone = parse("TSPO_traj.pdb")
    dssp = DSSP(backbone)
    list_helices = dssp.get_ca()
    for helix in list_helices:
        orig, axis = principal_axis(helix)
        len = length(helix, axis, orig)
        print(len)


test_length()
