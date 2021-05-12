def test_parse():
    # Test for static molecule
    backbone = parse("glut1.pdb")
    print("#a_carbon: ", len(backbone["alpha_carbon"]))
    print("#carbon:   ", len(backbone["carbon"]))
    print("#hydrogen: ", len(backbone["hydrogen"]))
    print("#nitrogen: ", len(backbone["nitrogen"]))
    print("#nitrogen: ", len(backbone["oxygen"]))
    print("#res:      ", len(backbone["res_number_list"]))

    # Test for trajectory file
    backbones = parse("TSPO_traj.pdb")
    backbone = backbones[2]
    print(len(backbones))
    print("#a_carbon: ", len(backbone["alpha_carbon"]))
    print("#carbon:   ", len(backbone["carbon"]))
    print("#hydrogen: ", len(backbone["hydrogen"]))
    print("#nitrogen: ", len(backbone["nitrogen"]))
    print("#nitrogen: ", len(backbone["oxygen"]))
    print("#res:      ", len(backbone["res_number_list"]))


# test_parse()
