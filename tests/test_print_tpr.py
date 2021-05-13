from context import tpr
from math import pi


def test_print_tpr_real():
    list_thetas = tpr.trajectory_tpr("data/TSPO_traj.pdb")
    list_tpr = [list_thetas[i][1] for i in range(len(list_thetas))]

    tpr.print_tpr(list_tpr)


def test_print_tpr_mock():
    list_tpr = [pi/2, 1.1*pi/2]
    tpr.print_tpr(list_tpr)


test_print_tpr_real()
