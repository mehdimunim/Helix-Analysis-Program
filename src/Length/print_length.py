from length_traj import trajectory_length
import matplotlib.pyplot as plt
import numpy as np


def print_length(list_lengths):
    """
    Print scatter plots of lengths

    """
    ind = [i for i in range(len(list_lengths))]

    list_lengths = np.array(list_lengths)
    ind = np.array(ind)

    plt.scatter(ind, list_lengths)
    plt.xlabel("frame numbers")
    plt.ylabel("length (A)")
    plt.title("Lengths of helix")
    plt.show()


def test_print_length_mock():
    list_lengths = [i for i in range(100)]
    print_length(list_lengths)


def test_print_length_real(filename):
    import os
    os.chdir(
        "C:\\Users\\Mehdi\\Documents\\GitHub\\Interdisciplinary-Project\\resource")
    list = trajectory_length("TSPO_traj.pdb")
    list_lengths = [list[i][1] for i in range(len(list))]
    print_length(list_lengths)


# test_print_length_mock()
test_print_length_real("TSPO_traj.pdb")
