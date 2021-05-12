def test_print_length_mock():
    list_lengths = [i for i in range(100)]
    print_length(list_lengths)


def test_print_length_real(filename):
    list = trajectory_length("TSPO_traj.pdb")
    list_lengths = [list[i][1] for i in range(len(list))]
    print_length(list_lengths)


test_print_length_real("TSPO_traj.pdb")
