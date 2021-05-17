def getAtomsPositions(filename):
    """ returns a dictionnary from the pdb file which contains coordinates of each atoms for the keys
        and its atom type for the value.
    """
    pdb = open(filename, "r")

    # déclaration du dictionnaire qui comportera toutes les informations
    res_coords = {}

    isBetween = False

    # Itérer sur les lignes dans pdb
    for i, line in enumerate(pdb):
        #  Vérifier si la ligne commence par "ATOM"
        if line.startswith('ATOM'):
            atom_name = line[12:16].strip()
            residueX = line[30:38].strip()
            residueY = line[38:46].strip()
            residueZ = line[46:54].strip()
            atom_name_short = line[76:100].strip()
            res_coords[(float(residueX), float(residueY),
                        float(residueZ))] = atom_name_short

    return res_coords
