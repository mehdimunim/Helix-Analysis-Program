import numpy as np

def dsmatvec(rot, din):
	"""
	Compute dout = din*rot (Vector*Matrix)
	***
	Parameters:
	rot: rotational matrix 3*3
	din: input vector 
	"""
	dout = np.zeros(3)
	for k in range(3):
		dout[k]=sum([din[i]*rot[k][i] for i in range(3)])
	return dout
