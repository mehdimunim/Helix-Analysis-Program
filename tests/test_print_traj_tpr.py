def test_traj_tpr(filename):
    import math
    list_thetas = trajectory_tpr(filename)
    print(list_thetas[0][1]*180/math.pi)
    print(list_thetas[1][1]*180/math.pi)
    print(list_thetas[2][1]*180/math.pi)
    print(list_thetas[3][1]*180/math.pi)
    print(list_thetas[4][1]*180/math.pi)
    print(list_thetas[5][1]*180/math.pi)


# test_traj_tpr("TSPO_traj.pdb")
