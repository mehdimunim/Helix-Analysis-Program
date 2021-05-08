def getAtomsPositions(filename):
    """ renvoie un dictionnaire à partir du fichier pdb qui contient les coordonnées de chaque atome
        en clé et son atome en valeur
    """
    pdb = open(filename, "r")

    #déclaration du dictionnaire qui comportera toutes les informations
    res_coords = {}

    # Itérer sur les lignes dans pdb
    for line in pdb:
        #  Vérifier si la ligne commence par "ATOM"
        if line.startswith('ATOM'):
            atom_name = line[12:16].strip()
            residueX = line[30:38].strip()
            residueY = line[38:46].strip()
            residueZ = line[46:54].strip()
            atom_name_short = line[70:100].strip()

            res_coords[(float(residueX),float(residueY),float(residueZ))] = atom_name_short

    return res_coords

def getAtomsPositionsBetween(filename,begin,end):
    """ renvoie un dictionnaire à partir du fichier pdb qui contient les coordonnées de chaque atome
        en clé et son atome en valeur. Ne renvoie que les atomes se trouvant entre begin et end
    """
    pdb = open(filename, "r")

    #déclaration du dictionnaire qui comportera toutes les informations
    res_coords = {}

    isBetween = False

    # Itérer sur les lignes dans pdb
    for line in pdb:
        #  Vérifier si la ligne commence par "ATOM"
        if line.startswith('ATOM'):
            atom_name = line[12:16].strip()
            residueX = line[30:38].strip()
            residueY = line[38:46].strip()
            residueZ = line[46:54].strip()
            atom_name_short = line[70:100].strip()
            if (float(residueX),float(residueY),float(residueZ)) == begin:
                isBetween=True
            if (float(residueX),float(residueY),float(residueZ)) == end:
                return res_coords
            if isBetween:
                res_coords[(float(residueX),float(residueY),float(residueZ))] = atom_name_short

    return res_coords

def getDifferentsAtoms(atomsCoords):
    """renvoie une liste de tous les différents atomes rencontrés dans le dictionnaire atomsCoords
       Juste pour voir les données que nous aurons besoin dans le dictionnaire 'atomicMass' 
    """
    atoms = set()
    for coo,atom in atomsCoords.items():
        atoms.add(atom)
    return atoms


def centreGravite(filename,begin,end):
    """Calcule le centre de gravité de l'hélice commençant aux coos begin et terminant aux coos end
       à partir de getAtomsPositionsBetween 
    """
    #On crée un dictionnaire pour relier un atome avec sa masse:
    atomsPositions = getAtomsPositionsBetween(filename,begin,end)

    # ATTENTION: les masses atomiques de O1- et N1+ ne sont peut être pas les bonnes
    atomicMass = {'C':12.0107,'H':1.00784,'O':15.999,'O1-':15.999,'S':16,'N':14.0067,'N1+':14.0067}

    xGrav,yGrav,zGrav = 0,0,0
    totalMass = 0

    for coo,name in atomsPositions.items():
        xGrav += coo[0]*atomicMass[name]
        yGrav += coo[1]*atomicMass[name]
        zGrav += coo[2]*atomicMass[name]
        totalMass+=atomicMass[name]
    
    return xGrav/totalMass,yGrav/totalMass,zGrav/totalMass





centreGravite("DSSP/glut1.pdb")