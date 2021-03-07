import math

def calcturnperres(nres,incrot,perpvec,axisdir,anglechangeref,irefang,MAXHX):
	"""

	Calculate TPR (Turn Per Residue)
	***
	Parameters:
	nres: number of residue
	incrot: increment rotation
	perpvec: perpendicular vector
	axisdir: axis vector
	anglechangeref: reference of angle change
	irefang; reference angle
	MAXHX: total number of helices

	Returns:
	turnperres : turn per residue
	
	"""
	PI = math.pi
	nresu = nres - 2*incrot
	a11 = (2*nresu**3 + 3*nresu**2 + nresu)/6
	a21 = (nresu*(nresu + 1))/2
	a12 = a21
	a22 = nresu
	b1 = 0
	b2 = 0
	turna = 0
	nturn = 0
	turnchangeav = 0
	for ir in range(1 + incrot, nres-incrot):
		turnchange = angcomp(perpvec[0][ir-1],axisdir,perpvec[0][ir])
		if (turnchange < 0):
			turnchange+=2*PI
		if (irefang == 1):
			#Compare with turn angle in the reference conformation
			if (abs(turnchange-anglechangeref(ir)) > abs((turnchange-2.0*PI)-anglechangeref(ir))):
				turnchange-= 2*PI
			print("ir = {} turnchange = {} ref= {} ".format(ir,turnchange*180/PI,anglechangeref[ir]*180/PI))
		turnchangeav+=turnchange
		turna+=turnchange
		print("ir={}, turna = {}, turnchange= {}".format(ir,turna*180/PI,turnchange*180/PI))
		b1+= (ir - incrot)*turna
		b2+=turna
	turnperres= (aa2*b1 - a12*b2)/(a22*a11 - a12*a21)
	turnchangeav=turnchangeav/(nres - 2*incrot -1)
	print("turnchangeav= {} {}",turnchangeav,turnchangeav*180/PI)
	return turnperres
	