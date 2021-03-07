import math
import numpy as np

def dacoscheck(arg):
	"""
	Make sure cos is in the [-1,1] range, return acos(arg)
	otherwise raise an exception
	***
	Parameters:
	dargin: input angle
	argin:
	idbp:

	Returns:
	acos(dargin) if idb=1 else acos(argin)
	"""
	#if (idb == 1):
		#arg = dargin
	#else:
		#arg = argin
	#Checking arg value and returning acos(arg) if correct
	if ((arg > 1 and arg<=1.01) or (arg < -1 and arg >=-1.01)):
		# arg is close enough to the [-1,1] range to be approximated by its boundaries
		arg = np.sign(arg)
	elif( arg >1.01 or arg < -1.01):
		raise ValueError("PROGRAM ERROR: invalid sin or cos value: " +arg)
	return math.acos(arg)