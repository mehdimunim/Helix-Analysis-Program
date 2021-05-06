
def alpha_carbons(filename):
    """
    Returns a list of tuples containing each:
    - the residue number
    - the coordinates of the alpha carbon
    """

    alpha_carbons = []
    res_number_list = []

    # Itérer sur les lignes dans pdbname
    with open(filename, "r") as pdbname:
        for line in pdbname:

            if line.startswith('ATOM'):

                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                res_number = line[23:26].strip()
                res_number = int(res_number)
                x = line[30:38].strip()
                y = line[38:46].strip()
                z = line[46:54].strip()

                atom_coordinates = (float(x), float(y), float(z))

                if res_name != "PRO" and atom_name == 'CA':
                    res_number_list.append(res_number)
                    alpha_carbons.append(atom_coordinates)

    list = []

    for nres, ca in zip(res_number_list, alpha_carbons):
        list.append((nres, ca))

    return list


def test_alpha_carbons():
    list = alpha_carbons("glut1.pdb")
    for tuple in list:
        print(tuple[0], tuple[1])


# test_alpha_carbons()
