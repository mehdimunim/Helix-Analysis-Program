import numpy as np
import math
import warnings

def check_rotmat(rot):
	"""
	Check if rot is a rotation matrix
	***
	Parameters:
	rot: rotation matrix
	"""
	# Need to make rotwarn, nrottest, devmax global???
	# global nrotwarn
	# global nrottest
	# global devmax
	u = np.array([[1,0,0], [0,1,0], [0,0,1]])
	# check the rotation matrix
	for i in range(3):
		for j in range(3):
			sum =0
			for k in range(3):
				sum+=rot[i][k]*rot[j][k]
			if (abs(sum - u[i][j]) >= 0.99):
				raise ValueError("Not a rotation matrix {} \n SUM rot({},k)*rot({},k)",rot,i,j)
			elif (abs(sum - u[i][j]) > 0.1):
				if(nrotwarn <= 0):
					warnings.warn("Not a rotation matrix {} \n SUM rot({},k)*rot({},k)",rot,i,j)
					nrotwarn+=1
				if (devmax < abs(sum - u[i][j])):
					devmax = abs(sum - u[i][j])
	print(" The rotation matrix in rotate_c:",rot)
   