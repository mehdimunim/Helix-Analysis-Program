from parse import parse
from length import length as get_length
from dssp import *
from axis import principal_axis


def trajectory_length(trajectory_file):
    """
    Get length list for each frame

    ---
    Parameters:
    trajectory_file: concatenation of pdb files

    ---
    Return:
    list_lengths: list of lists of lengths. 
    Each sub-list corresponds to a frame and each item of the sub-list to a helix

    """

    backbones = parse(trajectory_file)
    list_lengths = []
    for i, backbone in enumerate(backbones):
        dssp = DSSP(backbone)
        list_helices = dssp.get_ca()
        lengths = []
        for helix in list_helices:
            orig, axis = principal_axis(helix)
            length = get_length(helix, axis, orig)
            lengths.append(length)
        list_lengths.append(lengths)
    return list_lengths


def test_length_traj(filename):
    import math
    list_lengths = trajectory_length(filename)
    print(list_lengths[0][1])
    print(list_lengths[1][1])
    print(list_lengths[2][1])
    print(list_lengths[3][1])
    print(list_lengths[4][1])
    print(list_lengths[5][1])
    print(list_lengths[6][1])
    print(list_lengths[7][1])
    print(list_lengths[8][1])
    print(list_lengths[9][1])
    print(list_lengths[10][1])
    print(list_lengths[11][1])


# test_length_traj("TSPO_traj.pdb")
