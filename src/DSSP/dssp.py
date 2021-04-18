from parse import parse_structure
from hbonds import find_hbonds
from patterns import find_patterns


def dssp():
    """
    DSSP main function
    Inspired both from the paper of Kabsh & Sander and from Simulaid
    """
    pdb = "glut1.pdb"

    backbone = parse_structure(pdb)

    hbond = find_hbonds(backbone)

    secondary_structure = find_patterns(hbond)

    for type in [3, 4, 5]:

        name = str(type) + "-helices"

        print(" \n {}\n".format(name.upper()))

        str_starts = ""
        str_ends = ""
        for pos, val in enumerate(secondary_structure[name]):
            str_val = "{:4d}".format(val)
            if pos % 2 == 0:
                str_starts += str_val
            else:
                str_ends += str_val

        print("starts: ", str_starts)
        print("ends:   ", str_ends)

        print("Number of {}: {} ".format(
            name, int(len(secondary_structure[name])/2)))
