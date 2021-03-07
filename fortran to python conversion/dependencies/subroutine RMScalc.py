import math

def RMScalc(co, nats, dir, ip):
	"""

	Modified from the original to the calculate the sd of the atom to axis distance

	***
	Parameters
	co: the input coordinates
	nats: the number of atoms
	ip: the coordinates of the initial point
	dir: the axis vector

	Returns:
	the root mean squared of the set of point
	fp: the coordinates of the final point

	"""
	tmp = []
	# calculate rms deviations
	sum = 0
	sum2 = 0
	for i in range(nats):
		tmp = co[0][i] - ip
		tmp = dvproj(dir, tmp)
		fp = ip + tmp
		dd= ddistsq(co[0][i], fp)
		sum+= math.sqrt(dd)
		sum2+=dd
	return sqrt(sum2/nats - (sum/nats)**2),fp