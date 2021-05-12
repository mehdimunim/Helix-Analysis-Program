def test_hbonds():
    from parse import parse_structure
    backbone = parse_structure("glut1.pdb")

    hbonds = find_hbonds(backbone)
    print("#hbond", len(hbonds))
    print("example of hbond: ", hbonds[0])


# test_hbonds()
