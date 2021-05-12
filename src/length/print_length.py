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
