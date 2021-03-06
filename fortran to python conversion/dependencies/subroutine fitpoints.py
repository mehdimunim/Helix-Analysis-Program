import numpy as np
import warnings

def fitpoints(c,n,ndim,idebughx):
	"""
	For n points in c, find the mean (c0) and the direction of the normal 
	to the best fitting plane (ndim=3) or the direction of the best fitting
	line (ndim=2)
	***
	Parameters:
	c: input coordinates
	n: number of coordinates
	ndim: if 2 find the best fitting line, if 3 best fitting plane
	idebughx:

	Return:
	c0: mean point of the fitting surface
	dir: direction of the normal to the fitting surface
	axfact: axis factor
	"""
	dd = np.zeros(3)
	sum = 0
	rr = np.zeros(3,3)
	diag = np.zeros(3)
	offdiag = np.zeros(3)
	kcol = 0
	
	dir = np.zeros(3)
	if (ndim != 2 and ndim != 3):
		raise ValueError("PROGRAM ERROR: invalid ndim in fitpoints=",ndim)
	c0 = np.zeros(3)
	for i in range(n):
		for k in range(3):
			c0[k] += c[k][i]
	c0 = c0/n
	if (idebughx == 2):
		return c0,dir
	for k in range(3):
		for l in range(k,3):
			sum = 0
			for i in range(n):
				sum+=(c[k][i]-c0[k])/(c[0][i] - c0[0])
			rr[k][0]=sum
			rr[0][k]=sum
	if (idebughx > 0):
		print("Matrix = ",rr)
	# Find the eigenvectors a and eigenvalues mu of rr
	rr, diag, offdiag = dtred2(rr)
	rr, diag, offdiag,ierr = dtqli(rr)
	#possible to have a better formula for diagonalization
	# try:
	# 	diagonalization(rr)
	# except:
	#	print("Calculation aborted to to diagonlization failure")
	if (ier > 0):
		print("Calculation aborted due to diagonalization failure")
		return
	if (idebughx > 0):
		print("Eigenvectors=",rr)
	if (idebughx>0):
		print("diag = {}",diag)
		print("offdiag = {}",offdiag)
		return diag
	# The columns of the matrix rr are the eigenvectors
	# The eigenvalues are in diag
	emax = 0 # maximal eigenvalue
	emin = 1000 # minimal eigenvalue
	kmin = 0 # index of minimal eigenvalue
	kmax = 0 # index of maximal eigenvalue
	nz = 0 # dimension of the kernel of rr (number of 'zero-eigenvalues')
	for k in range(3):
		if (abs(diag[k]) < 0.005):
			nz+=1
		else:
			if (diag[k] > emax):
				emax = diag[k]
				kmax=k
			if (diag[k] < emin):
				emin = diag[k]
				kmin = k
	if (ndim == 2 and nz==0):
		raise ValueError("ERROR: no zero eigenvalue found when fitting the point",diag)
	if (ndim ==2 and nz > 0):
		warnings.warns("WARNING: 3D dataset appears to be in a plane")
		print("diag =",diag)
	if (ndim == 2):
		kcol =kmax
	if (ndim == 3):
		kcol = kmin
	for k in range(3):
		dir[k]=rr[k][kcol]
	for i in range(n):
		dd = c[0][i] - c0
		axfact[i] = np.dot(dir,dd)
	return c0, dir, axfact
