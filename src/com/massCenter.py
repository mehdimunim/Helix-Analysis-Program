from ..parser import getAtomsPositions
from ..parser import parse
from ..dssp import *
from ..parser import *
from numpy.linalg import *
import matplotlib.pyplot as plt
import numpy as np


def getDifferentsAtoms(atomsCoords):
    """renvoie une liste de tous les différents atomes rencontrés dans le dictionnaire atomsCoords
       Juste pour voir les données que nous aurons besoin dans le dictionnaire 'atomicMass'
    """
    atoms = set()
    for coo, atom in atomsCoords.items():
        atoms.add(atom)
    return atoms


def centreGravite(atomsPositions, helice):
    """Calcule le centre de gravité de l'hélice à partir de getAtomsPositionsBetween
    """
    # ATTENTION: les masses atomiques de O1- et N1+ ne sont peut être pas les bonnes
    atomicMass = {'C': 12.0107, 'H': 1.00784, 'O': 15.999,
                  'O1-': 15.999, 'S': 16, 'N': 14.0067, 'N1+': 14.0067}

    xGrav, yGrav, zGrav = 0, 0, 0
    totalMass = 0

    # to test if we are in the helice or not
    isBetween = False

    for coos, name in atomsPositions.items():

        if (coos == helice[0]):
            isBetween = True
        if (coos == helice[1]):
            return xGrav/totalMass, yGrav/totalMass, zGrav/totalMass

        if isBetween:

            xGrav += coos[0]*atomicMass[name]
            yGrav += coos[1]*atomicMass[name]
            zGrav += coos[2]*atomicMass[name]
            totalMass += atomicMass[name]

    return xGrav/totalMass, yGrav/totalMass, zGrav/totalMass


def getMassCenters(backbones, atomsPositions):
    res = []
    for backbone in backbones:
        dssp = DSSP(backbone)
        listHelices = dssp.get_ca()
        centresGravite = []
        # on suppose qu'une hélice est un tuple (begin=(x,y,z),end)
        for helice in listHelices:
            centresGravite.append(centreGravite(atomsPositions, helice))
        res.append(centresGravite)
    return res


def showGraphMassCenters_traj(filename):
    backbones = parse(filename)
    atomsPositions = getAtomsPositions(filename)
    res = getMassCenters(backbones, atomsPositions)

    # to get graphs:
    for massCenters in res:
        massCentersNormed = np.array([norm(x) for x in massCenters])
        axeX = [i for i in range(len(massCentersNormed))]
        plt.scatter(axeX, massCentersNormed)

    plt.show()


def showGraphMassCenters_static(filename):
    backbone = parse(filename)
    atomsPositions = getAtomsPositions(filename)
    res = getMassCenters([backbone], atomsPositions)

    # to get graphs:
    for massCenters in res:
        massCentersNormed = np.array([norm(x) for x in massCenters])
        axeX = [i for i in range(len(massCentersNormed))]
        plt.scatter(axeX, massCentersNormed)


def showGraphMassCenters(filename, isTrajectory):
    if isTrajectory:
        showGraphMassCenters_traj(filename)
    else:
        showGraphMassCenters_static(filename)
