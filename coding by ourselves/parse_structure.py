import numpy as np
import os
os.chdir(
    "C:\\Users\\Mehdi\\Documents\\GitHub\\Interdisciplinary-Project\\coding by ourselves")


def parse_structure(filename):
    """
    Parse PDB entry file
    to get the CA, C, O, N and H from the backbone chain
    Returns: the corresponding coordinates for each residue
    """
    alpha_carbons = []
    simple_carbons = []
    oxygens = []
    nitrogens = []
    hydrogens = []

    # Itérer sur les lignes dans pdbname
    # debug
    i = 0
    with open(filename, "r") as pdbname:
        for line in pdbname:

            #  Vérifier si la ligne commence par "ATOM"
            if line.startswith('ATOM'):

                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                x = line[30:38].strip()
                y = line[38:46].strip()
                z = line[46:54].strip()

                atom_coordinates = (float(x), float(y), float(z))

                if (res_name != "PRO"):
                    if atom_name == 'CA':
                        alpha_carbons.append(atom_coordinates)

                    elif atom_name == "C":
                        simple_carbons.append(atom_coordinates)

                    elif atom_name == "O":
                        oxygens.append(atom_coordinates)

                    elif atom_name == "N":
                        nitrogens.append(atom_coordinates)

                    elif atom_name == "HN":
                        hydrogens.append(atom_coordinates)

    # removing first residuecarbons, alpha-carbons and nitrogens
    # as there are no OH in it
    alpha_carbons = alpha_carbons[1:-1]
    simple_carbons = simple_carbons[1:-1]
    nitrogens = nitrogens[1:-1]

    # extra last oxygen
    oxygens = oxygens[1:]

    # extra first hydrogen
    hydrogens = hydrogens[:-1]

    return alpha_carbons, simple_carbons, oxygens, nitrogens, hydrogens


def test_parse():

    alpha_carbons, simple_carbons, oxygens, nitrogens, hydrogens = parse_structure(
        "glut1.pdb")
    print("alpha carbons", len(alpha_carbons))
    print("simple carbons", len(simple_carbons))
    print("oxygens", len(oxygens))
    print("nitrogens", len(nitrogens))
    print("hydrogens", len(hydrogens))


test_parse()
