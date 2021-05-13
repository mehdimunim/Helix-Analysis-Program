from context import length
from context import axis as axis_mod
from context import dssp as dssp_mod
from context import parser


def test_length():
    backbone = parser.parse("data/glut1.pdb")
    dssp = dssp_mod.DSSP(backbone)
    list_helices = dssp.get_ca()
    for helix in list_helices:
        orig, axis = axis_mod.principal_axis(helix)
        len = length.length(helix, axis, orig)
        print(len)


def test_length_traj():
    backbone = parser.parse("data/TSPO_traj.pdb")
    dssp = dssp_mod.DSSP(backbone)
    list_helices = dssp.get_ca()
    for helix in list_helices:
        orig, axis = axis_mod.principal_axis(helix)
        len = length.length(helix, axis, orig)
        print(len)


test_length()
