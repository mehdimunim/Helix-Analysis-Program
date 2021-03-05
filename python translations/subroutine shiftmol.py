def shiftmol(c, n, c0, fac):
	"""
	Center the molecule at c0
	***
	Parameters:
	c: molecule coordinates (dimension 3*n)
	c0: given point in space 
	"""
	print("SHIFTMOL n = {} fac = {} c0= {} ".format(,n,fac,c0))
	for i in range(n):
		for k in range(3):
			cnew[k][i] = c[k][i] + fac*c0[k]
	return cnew