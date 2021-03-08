def paramx(cent,ax,a):
	"""
	
	***
	Parameters:
	cent, ax: two 3D vectors
	a: factor

	"""
	for k in range(3):
		x[k] = cent[k] + a*ax[k]
	return x