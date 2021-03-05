import numpy as np

def angcomp(v1,v2,d0):
	"""
	Calculate the angle between v1 and v2
	Obtain the sign by requiring that d0 and d correspond to the z axis
	***
	Parameters:
	v1 and v2: two input vectors
	d0:
	"""
	arg = np.dot(v0,v)/dsqrt(np.dot(v0,v0)*np.dot(v,v))
	ang = dacoscheck(arg)
	v01 = np.cross(v0,v)
	if (np.dot(v01,d0) < 0):
		ang=-ang
	return ang	