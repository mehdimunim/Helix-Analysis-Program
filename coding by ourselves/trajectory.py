
import MDAnalysis as mda
u = mda.Universe("TSPO_traj.pdb") 

print(u.trajectory.n_frames)

for ts in u.trajectory: 
    print(u)

with mda.Writer("glut1.pdb", multiframe=True) as pdb:
    for ts in u.trajectory:
        pdb.write(u)
print(u)



"""install MDAnalysis"""
""" il faut install biopython=1.68 """