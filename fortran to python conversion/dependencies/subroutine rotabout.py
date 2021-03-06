import math
import warnings

def rotabout(v, v0, t, axis):
	""" 
	Rotate vector v about point v0 by t radians along the axis specified
	***
	Parameters:
	v: vector to be rotated
	v0: reference vector
	t: angle (in radians) for the rotation
	axis: axis specified in upper or lower case
	"""
	vnew = []
	st = math.sin(t)
	ct = math.cos(t)
	if(axis.lower() == 'x'):
		vnew[0] = v[0]
		vnew[1] = v0[1] + v[1]*ct - v0[1]*ct - v[2]*st + v0[2]*st
		vnew[2] = v0[2] + v[2]*ct - v0[2]*ct + v[1]*st - v0[1]*st
	elif (axis.lower() == 'y'):
		vnew[0] = v0[0] + v[1]*ct - v0[0]*ct - v[2]*st + v0[1]*st
		vnew[1] = v[1]
		vnew[2] = v0[2] + v[2]*ct - v0[2]*ct - v[0]*st + v0[0]*st 
	elif(axis.lower() == 'z'):
		vnew[0] = v0[0] + v[0]*ct - v0[0]*ct - v[1]*st + v0[1]*st
		vnew[1] = v0[1] + v[1]*ct - v0[1]*ct + v[0]*st - v0[0]*st
		vnew[2] = v[2]
	else:
		warnings.warn("No rotation because wrong axis specified")
		return v
	return vnew


