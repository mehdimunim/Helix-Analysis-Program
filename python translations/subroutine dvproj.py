import numpy as np

def dvproj(x, y):
	""" 
	compute z = the projection on x of y
	projxy requires that x be a unit vector
	formula projxy = (y dot x)*x
	"""
	x_magnitude = np.linalg.norm(x)
	# Check condition on magnitude
	if x_magnitude!=1:
		print("X is not a unit vector")
		print("Normalizing X just in case")
		x = x/x_magnitude
	return np.dot(y,x)*x