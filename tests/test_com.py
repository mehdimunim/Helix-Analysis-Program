from context import parser
from context import com
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm

# test de la fonction centreGravite:


def test_com():
    filename = "data/TSPO_traj.pdb"
    backbones = parser.parse(filename)
    atomsPositions = parser.getAtomsPositions(filename)
    res = com.getMassCenters(backbones, atomsPositions)
    print("res done !")

    # to get graphs:
    for massCenters in res:
        massCentersNormed = np.array([norm(x) for x in massCenters])
        axeX = [i for i in range(len(massCentersNormed))]
        plt.scatter(axeX, massCentersNormed)

    plt.show()


test_com()
