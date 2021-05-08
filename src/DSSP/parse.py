#
# Gets the backbone atoms
#

def parse_structure(filename):
    """
    Parse PDB entry file
    to get the CA, C, O, N and H from the backbone chain
    Returns: the corresponding coordinates for each residue
    ---
    NOTE: Supress first and last residues as they are electrically charged
    """
    alpha_carbons = []
    simple_carbons = []
    oxygens = []
    nitrogens = []
    hydrogens = []
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

                if (res_name != "PRO"):

                    if atom_name == 'CA':
                        res_number_list.append(res_number)
                        alpha_carbons.append(atom_coordinates)

                    elif atom_name == "C":
                        simple_carbons.append(atom_coordinates)

                    elif atom_name == "O":
                        oxygens.append(atom_coordinates)

                    elif atom_name == "N":
                        nitrogens.append(atom_coordinates)

                    elif atom_name == "HN":
                        hydrogens.append(atom_coordinates)

    # removing first residue carbons, alpha-carbons and nitrogens
    # as there are no OH in it
    alpha_carbons = alpha_carbons[1:-1]
    simple_carbons = simple_carbons[1:-1]
    nitrogens = nitrogens[1:-1]

    # update this to res_number_list
    res_number_list = res_number_list[1:-1]

    # removing first residue's extra oxygen
    oxygens = oxygens[1:]

    # removing last residue's extra hydrogen
    hydrogens = hydrogens[:-1]

    # Put it into a dictionary
    backbone = {
        "oxygen": oxygens,
        "carbon": simple_carbons,
        "alpha_carbon": alpha_carbons,
        "nitrogen": nitrogens,
        "hydrogen": hydrogens,
        "res_number_list": res_number_list
    }

    return backbone


def test_parse():
    backbone = parse_structure("glut1.pdb")
    print("#a_carbon: ", len(backbone["alpha_carbon"]))
    print("#carbon:   ", len(backbone["carbon"]))
    print("#hydrogen: ", len(backbone["hydrogen"]))
    print("#nitrogen: ", len(backbone["nitrogen"]))
    print("#nitrogen: ", len(backbone["oxygen"]))
    print("#res:      ", len(backbone["res_number_list"]))


# test_parse()
