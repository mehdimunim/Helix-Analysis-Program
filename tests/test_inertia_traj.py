from context import axis


def test_inertia_traj(filename):
    axis.inertia_traj(filename)


test_inertia_traj("data/TSPO_traj.pdb")
