from context import parser
from context import com
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm

# test de la fonction centreGravite:


def test_com():
    filename = "data/TSPO_traj.pdb"
    backbones = parser.parse(filename)
    com.showGraphMassCenters(filename, True)



test_com()
