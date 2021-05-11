#
# Gets the backbone atoms
#


def parse(filename):
    """
    Parse a pdb file to get all backbone atoms

    ---
    Parameters:
    filename: trajectory or static pdb file

    ---
    Returns:
    if trajectory, list of backbones for each frames
    else, backbone of the molecule

    """
    limits = get_frames_limits(filename)
    backbones = []
    for limit in limits:
        backbone = parse_structure(filename, limit[0], limit[1])
        backbones.append(backbone)
    if (len(backbones) == 1):
        return backbones[0]
    return backbones


def parse_structure(filename, start, end):
    """
    Parse PDB entry file
    to get the CA, C, O, N and H from the backbone chain

    ---
    Parameters:
    filename: pdb input file (trajectory or static)
    start: line from which to start reading the file
    end: line from which to stop reading

    ---
    Returns:
    backbone: the corresponding coordinates for each residue

    ---
    NOTE: Supress first and last residues
    """

    alpha_carbons = []
    simple_carbons = []
    oxygens = []
    nitrogens = []
    hydrogens = []
    res_number_list = []

    with open(filename, "r") as pdbname:
        # Select slice to read
        lines = pdbname.readlines()[start:end+1]
        for line in lines:

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

                    elif atom_name == "HN" or atom_name == "H":
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


def get_frames_limits(filename):
    """
    Get starts and ends of each frames a trajectory file
    Get start and end of the file for a static file
    ---
    Parameters:
    filename: trajectory or static pdb file

    ---
    Return:
    frames: list of tuple with (startline, endline) for each frame

    """
    frames = []
    with open(filename, "r") as pdbname:
        lines = pdbname.readlines()
        start = 0
        end = len(lines)
        for i, line in enumerate(lines):
            if line.startswith("MODEL"):
                start = i

            elif line.startswith("ENDMDL"):
                end = i
                frames.append((start, end))
    if len(frames) == 0:
        frames.append((start, end))
    return frames


def test_parse():
    import os
    os.chdir(
<<<<<<< Updated upstream
        "C:\\Users\\Mehdi\\Documents\\GitHub\\Interdisciplinary-Project\\resource")
=======
        "C:\\Users\\hassn\\Documents\\A 2A\\Github\\src")
>>>>>>> Stashed changes

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


#test_parse()
