from context import length


def test_print_length_mock():
    list_lengths = [i for i in range(100)]
    length.print_length(list_lengths)


def test_print_length_real(filename):
    list = length.trajectory_length(filename)
    list_lengths = [list[i][1] for i in range(len(list))]
    length.print_length(list_lengths)


test_print_length_real("data/TSPO_traj.pdb")
