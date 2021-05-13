from context import dssp as dssp_mod
from context import parser


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
