import math
import numpy as np

def dacoscheck(arg):
	"""
	Make sure cos is in the [-1,1] range, return acos(arg)
	otherwise raise an exception
	***
	Parameters:
	arg: input angle
	"""
	#Checking arg value and returning acos(arg) if correct
	if ((arg > 1 and arg<=1.01) or (arg < -1 and arg >=-1.01)):
		# arg is close enough to the [-1,1] range to be approximated by its boundaries
		arg = np.sign(arg)
	elif( arg >1.01 or arg < -1.01):
		raise ValueError("PROGRAM ERROR: invalid sin or cos value: " +arg)
	return math.acos(arg)