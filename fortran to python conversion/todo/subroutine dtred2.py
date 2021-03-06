import math
import numpy as np

def dtred2(a,n,np,d,e):
	"""
	Calculate ???
	***
	Parameters:
	a: dimension np*np
	n:??
	np:??
	d: diagonal np
	e: offdiagonal np
	"""
	scale = 0
	h = 0
	f = 0
	g =0
	if (n > 1):
		for i in range(n,1,-1):
			l=i-1
			h=0
			scale=0
			if(l > 1):
				for k in range(l):
					scale+=abs(a[i][k])
				if (scale == 0):
					e[i]=a[i][l]
				else:
					for k in range(l):
						a[i][k]=a[i][k]/scale
						h+=a[i][k]**2
					f = a[i][l]
					g = -np.sign(math.sqrt(h),f)
					e[i] = scale*g
					h-= f*g
					a[i][l]=f-g
					f=0
					for j in range(l):
						a[j][i]=a[i][j]/h
						g=0
						for k in range(j):
							g+=a[j][k]*a[i][k]
						if (l > j):
							for k in range(j+1,l):
								g+=a[k][j]*a[i][k]
						e[j]=g/h
						f+=e[j]*a[i][j]
					hh=f/2h
					for j in range(l):
						f = a[i][j]
						g = e[j]-hh*f
						e[j]=g
						for k in range(j):
							a[j][k]=a[j][k] - f*e[k]-g*a[i][k]
			else:
				e[i]=a[i][l]
			d[i]=h
	d[0]=0
	e[0]=0
	for i in range(n):
		l = i-1
		if (d[i] != 0):
			for j in range(l):
				g = 0
				for k in range(l):
					g+=a[i][k]*a[k][j]
				for k in range(l):
					a[k][j]-=g*a[k][i]
		d[i] = a[i][i]
		a[i][i]=1
		if (l > 0):
			for j in range(l):
				a[i][j]=0
				a[j][i]=0
	return a
