from parse import parse
from tpr import tpr
from dssp import *
from axis import principal_axis


def trajectory_tpr(trajectory_file):
    """
    Get tpr list for each frame

    ---
    Parameters:
    trajectory_file: concatenation of pdb files

    ---
    Return:
    thetas: list of list of tprs

    """

    backbones = parse(trajectory_file)
    thetas = []
    for backbone in backbones:
        list_helices = DSSP(backbone).get_ca()
        for helix in list_helices:
            orig, axis = principal_axis(helix)
            theta = tpr(helix, axis, orig)
            thetas.append(theta)
    return thetas


def test_traj_tpr(filename):
    import os
    os.chdir(
        "C:\\Users\\Mehdi\\Documents\\GitHub\\Interdisciplinary-Project\\resource")
    thetas = trajectory_tpr(filename)
    print(thetas[0][2])
    print(thetas[1][2])


test_traj_tpr("TSPO_traj.pdb")
