
from dssp import *
from axis2 import *

import numpy
import math
import MDAnalysis as mda

def test_angle():
    u = mda.Universe("glut1.pdb") 
    dssp = DSSP("glut1.pdb")
    for ts in u.trajectory: 
        list_helices = dssp.get_ca()
        helix_1 = list_helices[0]
        helix_2 = list_helices[1]
        axis1 = principal_axis2(helix_1)
        axis2 = principal_axis2(helix_2)
        dot_product = numpy.dot(axis1, axis2)
        angle = numpy.arccos(dot_product)
        print(math.degrees(angle))


test_angle()


