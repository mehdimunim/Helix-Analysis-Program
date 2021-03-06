import numpy as np

def circfit(x, nats, dir, ip):
	""" 
	Calculates the best-fit circle to a set of data and returns the 
	center in ip
	***
	Parameters
	x: the coordinates (array of 3,2*MAXHX length)
	nats: the number of atoms
	dir: input direction vector
	ip: output center of x

	"""
	# read input
	xp = np.zeros(3,100)
	xp0 = np.zeros(3)
	tmp = np.zeros(3)
	zero = np.zeros(3)
	
	print("got nats:",nats)
	for i in range(nats):
		for k in range(3):
			print("got coords:",x[k][i])

	r,thet,fi = polar(dir[1], dir[2], dir[3])
	print("theta,phi",thet,fi)
	for iat in range(nats):
		xp[0][iat] = x[0][iat].copy()
		xp[0][iat] = rotabout(xp[0][iat],zero, -thet,'z')
		xp[0][iat] = rotabout(xp[0][iat],zero,-fi,'y')
	tmp = dir.copy()
	sum = np.zeros(5)

	for i in range(nats):
		for j in range (i+1, nats):
			sum[0] += (xp[0][i] - xp[0][j])**2
			sum[1] += (xp[1][i] - xp[1][j])**2
			sum[2] += (xp[2][i] - xp[2][j])*(xp[1][i] -xp[1][j])
			scrat = xp[0][i]**2 + xp[1][i]**2 - xp[0][j]**2 - xp[1][j]**2
			sum[3] += scrat*(xp[0][i] - xp[0][j])
			sum[4] += scrat*(xp[1][i] - xp[1][j])
	print("sum=",sum)		
	xp0[0] = (sum[3]*sum[1] - sum[4]*sum[2])/(2*(sum[0]*sum[1] - sum[2]**2))
	xp0[1] = (sum[3] - 2*xp0[0]*sum[0])/(2*sum[2])
	xp0[2] = xp[2][0] #z coordinate of first atom
	xp0 = rotabout(xp0,zero,fi,'y')
	xp0 = rotabout(xp0,zero,thet,'z')
	ip = xp0.copy() 
	print("Output circle ip = {},dir = {}".format(ip,dir))
	return ip,dir
