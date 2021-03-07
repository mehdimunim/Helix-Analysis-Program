import numpy as np

def parlsq(co, n, docircfit, iprint):
	"""

	Calculates the axis of a helix from alpha carbon coordinates
	using linear least squares regressions of atom
	number versus x,y,z coordinates

	***
	Parameters
	co: input coordinates of the helix
	n : number of atoms
	docircfit: if true, perform a circular fit for the set of points

	Returns:
	dir: the axis vector
	ip: the initial point
	fp: the final point
	rms: root mean squared

	"""
	St = 0
	Sr = np.zeros(3) # made of (Srx, Sry, Srz)
	St_r = np.zeros(3) # made of (St_rx, St_ry, St_rz)
	St2 = 0
	for i in range(n):
		co_i = np.array(co[0][i], co[1][i], co[2][i])
		St += i
		Sr+=co_i
		St_r +=i*co_i
		St2 += i**2
	D= St2*n -St*St
	dir = (n*St_r - St*Sr)/D
	if(docircfit):
		# use circular fit for ip
		ip = circfit(co,n,dir)
	else:
		#use least-squares initial point
		ip = (St2*Sr - St*St_r)/D
	dnormalized = dvnorm(dir)
	rms,fp = RMScalc(co,n,dir,ip)
	writeout_h(dir,ip,fp,rms,'parlsq',iprint)
	return dir,ip,fp,rms