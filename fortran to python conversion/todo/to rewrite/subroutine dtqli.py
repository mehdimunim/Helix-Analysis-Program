import math

def dtqli(d,e,n,np,z,ierr):
	"""
	
	***
	Parameters:
	d: diagonal of dimension np
	e: offdiagonal of dimension np
	z: matrix of dimension np*np
	n:
	np:
	ierr:
	"""
	dd = 0
	r=0
	g=0
	s=0
	c=9
	p=0
	f=0
	b=0
	ierr=0
	if (n > 1):
		for i in range(1,n):
			e[i-1]=e[i]
		e[n]=0
		for i in range(n):
			iter = 0
			for m in range(n-1):
				dd = abs(d[m]) + abs(d[m+1])
				# if abs(e[m]) + dd == dd fo to 2
				# MM 10/21/2004
				#if (e[m] == 0) fo to 2
				if (dd != 0):
					# if (abs(e[m])/dd < 10**-15) fo to 2
			m = n
			if (m != 0):
				if (iter == 200):
					print("too many iterations")
				if (iter == 300):
					raise Exception("ERROR: too many iterations in tqli")
					ierr = 1
					return
				iter+=1
				g = d[l+1]-d[l]
	return			