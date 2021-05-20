from context import length


def test_length_traj(filename):
    length.length_traj_corrected(filename)


test_length_traj("data/TSPO_traj.pdb")  
