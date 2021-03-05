import numpy as np

def parlsq(co, n, docircfit,dir, ip,fp, rms, iprint):
	""" calculates the axis of a helix from alpha carbon coordinates
	using linear least squares regressions of atom
	number versus x,y,z coordinates
	***
	Parameters
	n : number of atoms
	dir : the axis vector
	ip : the coordinates of the initial point
	fp : the coordinates of the final point
	"""
	St = 0
	Sr = 0 # made of (Srx, Sry, Srz)
	St_r = 0 # made of (St_rx, St_ry, St_rz)
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
		circfit(co,n,dir,ip)
	else:
		#use least-squares initial point
		ip = (St2*Sr - St*St_r)/D
	dnormalized = dvnorm(dir)
	RMScalc(co,n,dir,ip,fp,RMS)
	#writeout_h(dir,ip,fp,rms,'parlsq',iprint)