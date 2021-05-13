from context import dssp as dssp_mod
from context import parser


def test_get_ca(dssp):
    cas = dssp.get_ca()
    print(cas[3])


def test_get_helices(dssp):
    all_helices = dssp.get_helices()
    print(all_helices[3])


def test_print(dssp):
    dssp.basic_print()


def test():
    backbone = parser.parse("data/glut1.pdb")
    dssp = dssp_mod.DSSP(backbone)
    test_print(dssp)
    #print("\n *********** GET CAS TEST *********** \n")
    # test_get_ca(dssp)
    print("\n *********** GET HELICES TEST *********** \n")
    test_get_helices(dssp)


test()
