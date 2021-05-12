from dssp import *
from parse import *
from numpy.linalg import *
from matplotlib import *
import numpy as np


def getAtomsPositionsBetween(filename, begin, end):
    """ renvoie un dictionnaire à partir du fichier pdb qui contient les coordonnées de chaque atome
        en clé et son atome en valeur. Ne renvoie que les atomes se trouvant entre begin et end
    """
    pdb = open(filename, "r")

    # déclaration du dictionnaire qui comportera toutes les informations
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
            atom_name_short = line[76:100].strip()
            if (float(residueX), float(residueY), float(residueZ)) == begin:
                isBetween = True
            if (float(residueX), float(residueY), float(residueZ)) == end:
                return res_coords
            if isBetween:
                res_coords[(float(residueX), float(residueY),
                            float(residueZ))] = atom_name_short

    return res_coords


def getDifferentsAtoms(atomsCoords):
    """renvoie une liste de tous les différents atomes rencontrés dans le dictionnaire atomsCoords
       Juste pour voir les données que nous aurons besoin dans le dictionnaire 'atomicMass'
    """
    atoms = set()
    for coo, atom in atomsCoords.items():
        atoms.add(atom)
    return atoms


def centreGravite(filename, helice):
    """Calcule le centre de gravité de l'hélice à partir de getAtomsPositionsBetween
    """
    # On crée un dictionnaire pour relier un atome avec sa masse:
    atomsPositions = getAtomsPositionsBetween(filename, helice[0], helice[-1])

    # ATTENTION: les masses atomiques de O1- et N1+ ne sont peut être pas les bonnes
    atomicMass = {'C': 12.0107, 'H': 1.00784, 'O': 15.999,
                  'O1-': 15.999, 'S': 16, 'N': 14.0067, 'N1+': 14.0067}

    xGrav, yGrav, zGrav = 0, 0, 0
    totalMass = 0

    for coo, name in atomsPositions.items():
        xGrav += coo[0]*atomicMass[name]
        yGrav += coo[1]*atomicMass[name]
        zGrav += coo[2]*atomicMass[name]
        totalMass += atomicMass[name]

    return xGrav/totalMass, yGrav/totalMass, zGrav/totalMass


def getMassCenters(backbones):
    res = []
    for backbone in backbones:
        dssp = DSSP(backbone)
        listHelices = dssp.get_ca()
        centresGravite = []
        # on suppose qu'une hélice est un tuple (begin,end)
        for helice in listHelices:
            centresGravite.append(centreGravite(filename, helice))
        res.append(centresGravite)
    return res


# test de la fonction centreGravite:
def COM(filename):
    res = getMassCenters(backbones)
    print("res done !")
    axeX = []
    axeY = []
    for i, massCenters in enumerate(res):
        axeX.append(i)
        massCentersNormed = [norm(x) for x in massCenters]
        axeY.append(massCentersNormed)
    axeX = np.array(axeX)
    axeY = np.array(axeY)
    plot.set_title("Center of Mass")
    plot.scatter(axeX, axeY)
    plot.savefig("Center of Mass.png")
