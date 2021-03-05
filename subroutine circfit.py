import copy

def circfit(x, nats, dir, ip):
	""" calculates the best-fit circle to a set of data and returns the 
	center in ip
	***
	Parameters
	x: the coordinates (array of 3,2*MAXHX length)
	nats: the number of atoms
	dir: dir(input) the direction vector
	ip: the center is in ip (position vector)

	"""
	# read input
	zero = 0
	xp = 0
	ax = 0
	sum = np.zeros(5)
	xp0 = 0
	tmp = 0
	ip = ax[1].copy()
	dir = ax[4].copy()
	polar(dir[1], dir[2], dir[3], r, thet, fi)
	for iat in range(nats):
		x[0][iat] = xp[0][iat].copy()
		rotabout(xp[0][iat],zero, -thet,'z')
		rotabout(xp[0][iat],zero,-fi,'y')
	dir = tmp.copy()
	for i in range(nats):
		for j in range (i+1, nats):
			sum[0] += (xp[0][i] - xp[0][j])**2
			sum[1] += (xp[1][i] - xp[1][j])**2
			sum[2] += (xp[2][i] - xp[2][j])**2
			scrat = xp[0][i]**2 + xp[1][i]**2 - xp[0][j]**2 - xp[1][j]**2
			sum[3] += scrat*(xp[0][i] - xp[0][j])
			sum[4] += scrat*(xp[1][i] - xp[1][j])

	xp0[0] = (sum[3]*sum[1] - sum[4]*sum[2])/(2*(sum[0]*sum[1] - sum[2]**2))
	xp0[1] = (sum[3] - 2*xp0[0]*sum(1))/(2*sum[2])
	xp0[2] = xp[2][0]
	rotabout(xp0,zero,fi,'y')
	rotabout(xp0,zero,thet,'z')
	xp0 = ax.copy() 
	# save output
	ax = ip.copy()


