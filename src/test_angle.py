
from dssp import *
from axis2 import *
from parse import parse
from length import length as get_length
from dssp import *
from axis import principal_axis
from parse import test_parse


import numpy
import math


def trajectory_angle(trajectory_file):
    """
    Gives the angle between the first two helices detected by dssp 

    ---
    Parameters:
    trajectory_file: concatenation of pdb files

    ---
    Return:
    angle  in degree 


    """

    backbones = parse(trajectory_file)
    list_lengths = []
    for i, backbone in enumerate(backbones):
        dssp = DSSP(backbone)
        list_helices = dssp.get_ca()
        angles = []
        helix_1 = list_helices[0]
        helix_2 = list_helices[1]
        axis1 = principal_axis2(helix_1)
        axis2 = principal_axis2(helix_2)
        dot_product = numpy.dot(axis1, axis2)
        angle = numpy.arccos(dot_product)
    return math.degrees(angle)



def test_angle(filename):
    import os
    import math
    os.chdir(
        "C:\\Users\\hassn\\Documents\\A 2A\\Github\\src")
    angle = trajectory_angle(filename)
    print("Angle between the two first helices :")
    print(angle)


test_angle("TSPO_traj.pdb")
