from parse_structure import parse_structure
from H_bonds import find_Hbonds
from find_patterns import find_patterns
from check_irregularities import check_irregularities


def main():
    """
    DSSP main function
    """

    pdb_cleaned = check_irregularities("glut1")

    alpha_carbons, simple_carbons, oxygens, nitrogens, hydrogens = parse_structure(
        pdb_cleaned)

    hbond = find_Hbonds(
        alpha_carbons, simple_carbons, oxygens, nitrogens, hydrogens)

    secondary_structure = find_patterns(hbond)

    for type in [3, 4, 5]:

        name = str(type) + "-helices"

        print(" \n {} [start, end, start, end...]\n".format(name.upper()))

        print(secondary_structure[name])

        print("Number of {} {} ".format(
            name, int(len(secondary_structure[name])/2)))


main()
exit(0)
