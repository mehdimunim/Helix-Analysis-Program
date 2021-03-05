import numpy as np

def fitpoints(c,n,ndim,c0,dir,axfact,idebughx):
	"""
	For n points in c, find the mean (c0) and the direction of the normal 
	to the best fitting plane (ndim=3) or the direction of the best fitting
	line (ndim=2)
	"""
	dd = np.zeros(3)
	sum = 0
	ddot = 0
	rr = np.zeros(3,3)
	diag = np.zeros(3)
	offdiag = np.zeros(3)
	kcol = []
	kmin = []
	kmax = []
	if (ndim != 2 and ndim != 3):
		print("PROGRAM ERROR: invalid ndim in fitpoints= {}",ndim)
		return
	c0 = np.zeros(3)
	for i in range(n):
		for k in range(3):
			c0[k] += c[k][i]
	for k in range(3):
		c0[k]=c0[k]/n
	if (idebughx == 2):
		return c0
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
	dtred2(rr,3,3,diag,offdiag)
	dtqli(diag,offdiag,3,3,rr,ierr)
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
	emax = 0
	emin = 1000
	nz = 0
	for k in range(3):
		if (dabs(diag[k]) < 0.005):
			nz+=1
		else:
			if (diag[k] > emax):
				emax = diag[k]
				kmax=k
			if (diag[k] < emin):
				emin = diag[k]
				kmnin = k
	if (ndim == 2 and nz==0):
		print("ERROR: no zero eigenvalue found when fitting the point",diag)
		return diag
	if (ndim ==2 and nz > 0):
		print("WARNING: 3D dataset appears to be in a plane")
	if (ndim == 2):
		kcol =kmax
	if (ndim == 3):
		kcol = kmin
	for k in range(3):
		dir[k]=rr[k][kcol]
	for i in range(n):
		dd = c[0][i] - c0
		axfact[i] = np.dot(dir,dd)
	return
