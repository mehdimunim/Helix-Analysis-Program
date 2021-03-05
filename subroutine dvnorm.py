import numpy as np

def dvnorm(x):
	""" 
	normalize the vector x
	"""
	x_magnitude = np.linalg.norm(x)
	# check if x_magnitude is really close to zero
	if x_magnitude< 10**-12: 
       print("DVNORM: Vector x is null")
       return x
    return x / x_magnitude