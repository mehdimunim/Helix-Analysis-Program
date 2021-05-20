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


def length_list_corrected(trajectory_file):
    """
    Correction: apply dssp's first assignement to all frames
    Get list of helices lengths for each frame.

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
    # DSSP for first frame
    dssp_init = dssp_mod.DSSP(backbones[0])
    first_assignement = dssp_init.get_structures()
    for i, backbone in enumerate(backbones):
        dssp = dssp_mod.DSSP(backbone)
        # get CAs of current frame according to first frame assignement
        list_helices = dssp.get_ca_with(first_assignement)
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
    nhelices = len(list_helices)
    list_lengths = []
    for helix in list_helices:
        orig, axis = principal_axis(helix)
        len_ = length_helix(helix, axis, orig)
        list_lengths.append(len_)

    textstr = "Number of helices " + str(nhelices)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.text(0.05, 0.05, textstr, fontsize=12,
             verticalalignment='bottom', bbox=props)
    plt.xlabel("length (A)")
    plt.ylabel("Number of helices")
    plt.title("Distribution of helices' lengths in " + molecule_name)
    plt.hist(
        list_lengths, bins=nhelices, facecolor='red', alpha=0.5)

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


def length_traj_corrected(filename):
    """
    Correction: use dssp's first assignement
    Save scatter plots of lengths in output/length_helix_NUM_MOLECULE.png where NUM is helix number
    and MOLECULE the name of the molecule

    """
    # get molecule name from corrected file
    molecule_name = filename.split("/")[1][:-4]

    # get list of lengths for each frame and each helix
    list_lengths = length_list_corrected(filename)

    nhelix = min([len(frame) for frame in list_lengths])

    nframe = len(list_lengths)

    plt.close()

    for i in range(1, nhelix):

        traj = [frame[i-1] for frame in list_lengths]

        textstr = "Total number of frames " + str(nframe)
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        plt.text(0.05, 0.05, textstr, fontsize=12,
                 verticalalignment='bottom', bbox=props)
        plt.xlabel("length (A)")
        plt.ylabel("Number of frames")

        plt.title("Distribution of lengths for " +
                  "helix " + str(i) + " in " + molecule_name)

        plt.hist(
            traj, bins=nframe, facecolor='red', alpha=0.5)

        plt.savefig("output/length_helix_" + str(i) +
                    "_" + molecule_name + ".png")
        plt.clf()
