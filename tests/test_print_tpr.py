def test_print_tpr_real():
    list_thetas = trajectory_tpr("TSPO_traj.pdb")
    list_tpr = [list_thetas[i][1] for i in range(len(list_thetas))]

    print_tpr(list_tpr)


def test_print_tpr_mock():
    import math
    list_tpr = [math.pi/2, 1.1*math.pi/2]
    print_tpr(list_tpr)


# test_print_tpr_real()
