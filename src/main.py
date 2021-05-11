<<<<<<< Updated upstream
if __name__ == "main":
    print("Not yet implemented")
=======
from parse import parse
from length import length as get_length
from dssp import *
from axis import principal_axis
from parse import test_parse
from length_traj import test_length_traj
from trajectory_tpr import test_traj_tpr
from test_angle import test_angle


if __name__ == "main":
 
   test_traj_tpr("TSPO_traj.pdb")
   
   test_length_traj("TSPO_traj.pdb")
   test_angle("TSPO_traj.pdb")
   
>>>>>>> Stashed changes
