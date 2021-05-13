from context import tpr
import math


def test_traj_tpr(filename):
    list_thetas = tpr.trajectory_tpr(filename)
    print(list_thetas[0][1]*180/math.pi)
    print(list_thetas[1][1]*180/math.pi)
    print(list_thetas[2][1]*180/math.pi)
    print(list_thetas[3][1]*180/math.pi)
    print(list_thetas[4][1]*180/math.pi)
    print(list_thetas[5][1]*180/math.pi)


test_traj_tpr("data/TSPO_traj.pdb")
