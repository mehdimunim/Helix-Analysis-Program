from context import tpr
import math


def test_traj_tpr(filename):
    tpr.tpr_traj_corrected(filename)


test_traj_tpr("data/TSPO_traj.pdb")
