import numpy as np
import math

def checkbend (c,axisdir,camod,perpvec,n,axtol,idebughx):
	"""
	Check if strucutre c is bend
	***
	Parameters:
	c: input coordinates of dimension 3*n
	axisdir: axis vector
	camod: of dimension 3*n
	perpvec: perpendicular vector of dimension 3*n
	axtol: axis tolerance
	idebughx:

	Returns:
	idecide: integer corresponding to the type of bend: too short, bent, random, alternate, too long 
	nup,ndown,nrun,nnear: 
	x0,rcirc: center and radius of the structure
	rn : direction vector
	"""
	rm = np.zeros(3)
	x = np.zeros(3)
	x1 = np.zeros(3)
	x2 = np.zeros(3)
	rr = 0
	rmin = 0
	ddistsq = 0
	x0, rn, axfact = fitpoints(c,n,3,x0,rn,idebughx)
	
	if (idebughx > 0):
		for i in range(n):
			writepdbd(77,c[0][i],i,i,'C','CA','A',1,0)
		writepdbd(77,x0,n+1,n+1,'S','X0','B',1.0,0.0)
		if (idebughx > 1):
			idecide = 6
			nup = 0
			ndown = 0
			nrun = 0
		else:
			x1 = x0 + rn
			writepdbd(77,x1,n+1,n+2,"S","RN","B",1,0)
	
	if (idebughx > 0):
		for i in range(n):
			for k in range(3):
				print("{} c= {} camod={} axfac={} perp= {}".format(i,c[k][i], camod[k][i], axfact[i], perpvec[k][i]))
				print("\n")
	# Find the shortest axis-Calpha distance
	rmin = 1000
	for i in range(n):
		rr = ddistsq(c[0][i],camod[0][i])
		if (rr < rmin):
			rmin =rr
			imin = i
	rr = dsqrt(rmin)
	# "Pull" uniformly all the Calphas toward the axis as close as possible
	for i in range(n):
		for k in range(3):
			camod[k][i] = c[k][i] + rr*perpvec[k][i]
	
	if (idebughx > 0):
		for i in range(n):
			writepdbd(77,camod[0][i],i,i,'C','CA ','B',1.0,0.0)
	# Calculate the radius of the fitting circle to the pulled points
	circ = circfit(c,n)
	rcirc = math.sqrt(ddistsq(circ,x0))
	rr = 0
	for i in range(n):
		rr+= ddistsq(circ,c[0][i])
	rcirc = math.sqrt(rr/n)
	
	if(idebughx > 0):
		print("CHECKBEND circ= {} x0={} rcirc={}".format(circ,x0,rcirc))
	# Find the projection of the pulled points on the fitting plane
	for i in range(n):
		x1 = camod[0][i] - x0
		rr = np.dot(x1,rn)
		for k in range(3):
			camod[k][i]-= rr*rn[k]
	x0,rm,axfact = fitpoints(camod,n,2,0)
	
	if (idebughx > 0):
		for i in range(n):
			writepdbd(77,camod(1,i),i,i,'N','CAP','C',1.0,0.0)
		writepdbd(77,x0,n+1,n+1,'P','X0','B',1.0,0.0)
		x1 = x0 + rm
		writepdbd(77,x1,n+2,n+2,'P','RM','B',1.0,0.0)
	
	if (idebughx > 0):
		#Check if the projections are indeed in a plane
		x1 = camod[0][0] - x0
		x2 = camod[0][1] - x0
		x = dvprd(x1,x2)
		for i in range(2,n):
			x1 = camod[0][i] - x0
			rr = np.dot(x,x1)
			for k in range(3):
				print("{} camod={} plane check(->0)={}".format(i,camod[k][i],rr))
				print("\n")
		cx = camod.copy()

	# Get the distance vectors between the point and its projection on the line
	for i in range(n):
		for k in range(3):
			camod[k][i] = camod[k][i]-x0[k] - axfact[i]*rm[k]
	
	if (idebughx > 0):
		for i in range(n):
			x = cx[0][i] - camod[0][i]
			writepdbd(77,x,i,i,'O','CA0','D',1.0,0.0)
	nnear = 0
	if (axtol > 0):
		#Eliminate (and count) residuea that are within tolerance of the axis
		for i in range(n):
			if (np.dot(camod[0][i],camod[0][i]) <= axtol**2):
				nnear+=1
			elif (nnear > 0):
				camod[0][i-nnear] = camod[0][i].copy()
	nupdown[0]=0
	nupdown[1]=0
	iupdown = 0
	nswitch = 0
	for i in range(1,n -nnear):
		if (np.dot(camod[0][i-1],camod[0][i]) < 0):
			nswitch+=1
			iupdown= 1-iupdown
		nupdown[iupdown +1]= nupdown[iupdown+1]+1
	nrun = nswitch +1
	nup= nupdown[0]
	ndown = nupdown[1]
	idecide = runtest(nup,ndown,nrun,idecide)
	return idecide,nup,ndown,nrun,nnear, x0, rcirc, rn
