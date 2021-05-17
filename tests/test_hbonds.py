from context import dssp as dssp_mod
from context import parser

def test_hbonds():
    backbone = parser.parse("data/glut1.pdb")

    hbonds = dssp_mod.find_hbonds(backbone)
    print("#hbond", len(hbonds))
    print("example of hbond: ", hbonds[0])


test_hbonds()
