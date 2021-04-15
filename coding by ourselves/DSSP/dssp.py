from parse import parse_structure
from hbonds import find_hbonds
from patterns import find_patterns
from irreg import check_irregularities


def main():
    """
    DSSP main function
    """

    pdb_cleaned = check_irregularities("glut1")

    alpha_carbons, simple_carbons, oxygens, nitrogens, hydrogens = parse_structure(
        pdb_cleaned)

    hbond = find_hbonds(
        alpha_carbons, simple_carbons, oxygens, nitrogens, hydrogens)

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


main()
exit(0)
