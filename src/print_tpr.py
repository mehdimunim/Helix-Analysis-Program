import numpy as np
import matplotlib.pyplot as plt
from trajectory_tpr import trajectory_tpr


def print_tpr(filename):
    """
    Print turn angle per residue for all frames

   """
   list_thetas = trajectory_tpr(filename)
   len =  min([len(list) for list in list_thetas])
   for j in range(len(list_thetas[0])):
       list_tpr = [list_thetas[i][j] for i in range(len(list_thetas))]
       print_tpr(list_tpr, j)




def print_helix_tpr(list_tpr, nhelix):
    """
    Print turn angle per residue of the given frame as a dial plot

    ---
    Parameters:
    list_thetas: output of trajectory tpr (in rad)

    ---
    Output:
    Dial plot of the tpr for a given helix

    """

    plt.axes

    r = np.arange(0, 1, 1/100)

    fig, ax=plt.subplots(subplot_kw={'projection': 'polar'})

    ax.set_rticks([0.5, 0.75, 1])

    for theta in list_tpr:
        theta_to_plot = [theta for _ in r]
        ax.plot(theta_to_plot, r)

    ax.grid(True)

    ax.set_title("Turn angle per residue " + str(nhelix), va='bottom')
    plt.show()
    plt.savefig("Turn angle per residue " + str(nhelix) + ".png")


def test_print_tpr_real():
    list_thetas = trajectory_tpr("TSPO_traj.pdb")
    list_tpr = [list_thetas[i][1] for i in range(len(list_thetas))]

    print_tpr(list_tpr)


def test_print_tpr_mock():
    import math
    list_tpr = [math.pi/2, 1.1*math.pi/2]
    print_tpr(list_tpr)


# test_print_tpr_real()
