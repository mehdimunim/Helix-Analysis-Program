import numpy as np

def project(x,c,ax):
	"""
	
	Find the nearest point p from x on c+a*ax
	
	***
	Parameters:
	x, c, ax: three vectors

	Returns:
	p: the nearest point 
	a : ax|(x-c)

	"""
	a = np.dot(ax,x-c)
	for k in range(3):
		p[k] = c[k] + a*ax[k]
	return p,a 