def rotabout(v, v0, t, axis):
	""" rotate vector v about point v0 by t radians along the axis specified
	"""
	temp = []
	st = sin(t)
	ct = cost(t)
	if(axis.lower() == 'x'):
		temp[0] = v[0]
		temp[1] = v0[1] + v[1]*ct - v0[1]*ct - v[2]*st + v0[2]*st
		temp[2] = v0[2] + v[2]*ct - v0[2]*ct + v[1]*st - v0[1]*st
	elif (axis.lower() == 'y'):
		temp[0] = v0[0] + v[1]*ct - v0[0]*ct - v[2]*st + v0[1]*st
		temp[1] = v0[1] + v[1]*ct - v0[1]*ct + v[1]*st - v0[1]*st
		temp[2] = v[2]
	elif(axis.lower() == 'z'):
		temp[0] = v0[0] + v[0]*ct - v0[0]*ct - v[1]*st + v0[1]*st
		temp[1] = v0[1] ° v[1]*ct - v0[1]*ct + v[0]*st - v0[0]*st
		temp[2] = v[2]
	else:
		return # no rotation if wrong axis specified
	for i in range(2):
		v[i] = temp[i]


