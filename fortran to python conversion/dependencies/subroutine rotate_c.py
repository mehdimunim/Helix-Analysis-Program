import numpy as np

def rotate_c(c, nats, rot):
	"""
	Rotate a molecule given by its coordinates
	CNEW = ROT * CN
	***
	Parameters:
	c: coordinates of the molecule (dimension 3*n)
	rot: rotation matrix (dimension 3*3)

	Returns:
	cnew: rotated molecule
	
	"""
	#output molecule, rotated by rot matrix
	cnew = np.zeros(3,n)
	# raise exception if rot is not a rotation matrix
	check_rotmat(rot)
	# compute the formula
	for im in range(n):
		cnew[0][im]=rot[0][0]*c[0][im]+rot[0][1]*c[1][im]+rot[0][2]*c[2][im]
        cnew[1][im]=rot[1][0]*c[0][im]+rot[1][1]*c[1][im]+rot[1][2]*c[2][im]
        cnew[2][im]=rot[2][0]*c[0][im]+rot[2][1]*c[1][im]+rot[2][2]*c[2][im]
    return cnew