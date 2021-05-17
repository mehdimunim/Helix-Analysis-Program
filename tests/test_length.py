from context import length as length_mod
from context import axis as axis_mod
from context import dssp as dssp_mod
from context import parser


def test_length():
    backbone = parser.parse("data/glut1.pdb")
    dssp = dssp_mod.DSSP(backbone)
    list_helices = dssp.get_ca()
    for helix in list_helices:
        orig, axis = axis_mod.principal_axis(helix)
        len = length_mod.length_helix(helix, axis, orig)
        print(len)


test_length()
