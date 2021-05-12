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



