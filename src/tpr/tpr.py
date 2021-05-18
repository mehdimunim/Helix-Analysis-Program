#
# Turn per residue
#

from ..parser import parse
from ..axis import principal_axis
from ..dssp import dssp as dssp_mod
import numpy as np
import matplotlib.pyplot as plt
from ..common import normal
from ..common import angle
from .fit import fit


def tpr_algo(alpha_carbons, axis_direction, axis_center):
    """
    Calculate the turn angle per residue
    Inspired from TRAJELIX from simulaid.

    Two steps:

    1- Calculate the turn angles phi_i between first residue and i

    2- Fit phi_i = theta*i + phi_0 where theta is the turn angle per residue

    (Mezei, Filizola, 2006)

    ---
    Parameters:
    alpha_carbons: alpha carbons (used to calculate the angle)
    axis_direction: inertia axis direction of the structure (often an alpha-helix)
    axis_center: center of the intertia axis

    ---
    Return:
    theta : turn per angle per residue (in rad)

    """

    # First step
    # Calculate the turn angles between first residue and i

    phi = []
    phi_i = 0
    for i in range(1, len(alpha_carbons)):
        # calculate normal vectors for residue i-1 and i
        vec_before, _ = normal(axis_direction, axis_center, alpha_carbons[i-1])
        vec_after, _ = normal(axis_direction, axis_center, alpha_carbons[i])

        # angle between i-1 and i
        angle_between = angle(vec_before, vec_after)

        # angle between first residue and i
        phi_i += angle_between

        phi.append(phi_i)

    # Second step
    # Find turn angle per residue with linear regression

    theta = fit(phi)

    return theta


def tpr_helix(list_tpr, nhelix):
    """
    Print turn angle per residue of the given frame as a dial plot
    Save the result to PNG file

    ---
    Parameters:
    list_thetas: output of trajectory tpr (in rad)

    ---
    Output:
    Dial plot of the tpr for a given helix

    """

    plt.axes

    r = np.arange(0, 1, 1/100)

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})

    ax.set_rticks([0.5, 0.75, 1])

    for theta in list_tpr:
        theta_to_plot = [theta for _ in r]
        ax.plot(theta_to_plot, r)

    ax.grid(True)

    ax.set_title("Turn angle per residue " + str(nhelix), va='bottom')
    plt.savefig("output/TPR_" + str(nhelix) + ".png")


def tpr_trajectory(filename):
    """"
    Calculate tpr variations for each helix
    Save all dial graphs in ouput
    """

    list_tprs = tpr_list(filename)

    nhelix = min([len(frame) for frame in list_tprs])

    for i in range(nhelix):
        tpr_helix(list_tprs, i)


def tpr(list_helices, filename):
    """"
    Calculate tpr for helices in backbone
    """
    plt.axes

    molecule_name = filename.split("/")[1][:-4]

    thetas = []

    for helix in list_helices:
        orig, axis = principal_axis(helix)
        theta = tpr_algo(helix, axis, orig)
        thetas.append(theta)

    r = np.arange(0, 1, 1/100)
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.set_rticks([0.5, 0.75, 1])

    for theta in thetas:
        theta_to_plot = [theta for _ in r]
        ax.plot(theta_to_plot, r)

    ax.grid(True)

    ax.set_title("Turn angle per residue " + molecule_name, va='bottom')
    plt.savefig("output/TPR_" + molecule_name + ".png")


def tpr_list(trajectory_file):
    """
    Get tpr list for each frame

    ---
    Parameters:
    trajectory_file: concatenation of pdb files

    ---
    Return:
    list_thetas: list of lists of tprs. 
    Each sub-list corresponds to a frame and each item of the sub-list to a helix

    """

    backbones = parse(trajectory_file, False)
    list_thetas = []
    for i, backbone in enumerate(backbones):
        dssp = dssp_mod.DSSP(backbone)
        list_helices = dssp.get_ca()
        thetas = []
        for helix in list_helices:
            orig, axis = principal_axis(helix)
            theta = tpr_algo(helix, axis, orig)
            thetas.append(theta)
        list_thetas.append(thetas)
    return list_thetas
