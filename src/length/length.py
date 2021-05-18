from ..common import normal
from ..common import distance
from ..axis import principal_axis
from ..parser import parse
from .. import dssp as dssp_mod
import matplotlib.pyplot as plt
import numpy as np


def length_helix(helix, axis, orig):
    """
    Calculate the length of the helix

    Can serve to learn about compression and extension of a helix during MD simulation
    cf. TRAJELIX 
    """
    _, first_orig = normal(axis, orig, helix[0])
    _, last_orig = normal(axis, orig, helix[-1])

    return distance(first_orig, last_orig)


def length_list(trajectory_file):
    """
    Get list of helices lengths for each frame

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
        dssp = dssp_mod.DSSP(backbone)
        list_helices = dssp.get_ca()
        lengths = []
        for helix in list_helices:
            orig, axis = principal_axis(helix)
            length = length_helix(helix, axis, orig)
            lengths.append(length)
        list_lengths.append(lengths)
    return list_lengths


def length(list_helices, filename):
    """
    For static PDB file
    Print scatter plots of lengths of each helices in molecule

    ---
    Output:
    Save result in output/length_MOLECULE.png

    """
    molecule_name = filename.split("/")[1][:-4]
    list_lengths = []
    for helix in list_helices:
        orig, axis = principal_axis(helix)
        len_ = length_helix(helix, axis, orig)
        list_lengths.append(len_)

    # number of helices
    nhelix = len(list_lengths)
    ind = [i for i in range(nhelix)]

    list_lengths = np.array(list_lengths)
    ind = np.array(ind)

    plt.scatter(ind, list_lengths)
    plt.xlabel("#Helix")
    plt.ylabel("length (A)")
    plt.title("Lengths of helices in molecule")
    plt.savefig("output/length_" + molecule_name + ".png")


def length_traj(filename):
    """
    Save scatter plots of lengths in output/MOLECULE_length_helix_NUM.png where NUM is helix number

    """
    list_lengths = length_list(filename)

    nhelix = min([len(frame) for frame in list_lengths])

    nframe = len(list_lengths)

    for helix in range(nhelix):
        traj = []
        for frame in range(nframe):
            traj.append(list_lengths[frame][helix])
        ind = [i for i in range(nframe)]
        list_lengths = np.array(list_lengths)
        ind = np.array(ind)
        plt.scatter(ind, traj)
        plt.xlabel("frame numbers")
        plt.ylabel("length (A)")
        plt.title("Lengths of helix " + str(helix))
        plt.savefig("output/length_helix_" + str(helix) + ".png")
