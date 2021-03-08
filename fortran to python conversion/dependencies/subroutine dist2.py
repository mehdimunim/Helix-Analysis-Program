def dist2(c1,c2):
	"""
	Calculate distance between c1 and c2
	***
	Parameters:
	c1, c2: two 3D vectors 
	"""
	dist2 = 0
	for k in range(3):
		dist2 += (c1[k] -c2[k])**2
	return dist2
