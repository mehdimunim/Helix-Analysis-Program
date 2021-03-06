import math

def polar(x, y, z):
	"""
	Return polar coordinates for point (x, y, z) 
	NOTE: Unlike most conventions, phi is the angle between R and the Z axis
	***
	Parameters:
	x,y,z: point
	"""
	PI = math.pi
	r = (x**2 + y**2 + z**2)**0.5
	if (x!=0 and y!= 0):
		theta = math.atan(abs(y/x))
		if ( y > 0):
			if (x < 0):
				theta = PI - theta # quadrant 2
		else:
			if (x < 0):
				theta = PI + theta # quadrant 3
			if (x > 0):
				theta = 2*PI-theta # quadrant 4
	else:
		if (x == 0):
			if (y >= 0):
				theta = PI/2
			else:
				theta = 3*PI/2
		if (y == 0):
			if (x >= 0):
				theta = 0
			else:
				theta = PI
	if ( r == 0):
		phi = 0
	else:
		cosa = z/r
		phi = dacoscheck(cosa, ccc, 1,6, 'POLAR')
	return r,theta,phi

