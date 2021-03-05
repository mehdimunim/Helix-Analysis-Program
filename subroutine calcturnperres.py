import math

def calcturnperres(turnperres,nres,incrot,perpvec,axisdir,anglechangeref,irefang,MAXHX):
	nresu = nres - 2*incrot
	pi = math.pi
	a11 = (2*nresu**3 + 3*nresu**2 + nresu)/6
	a21 = (nresu*(nresu + 1))/2
	a12 = a21
	a22 = nresu
	b1 = 0
	b2 = 0
	turna = 0
	nturn = 0
	turnchangeav = 0
	for ir in range(2 + incrot, nres-incrot):
		turnchange = angcomp(perpvec[0][ir],axisdir,perpvec[0][ir])
		if (turnchange < 0):
			turnchange+=2*pi
		if (irefang == 1):
			#Compare with turn angle in the reference conformation
			if (abs(turnchange-anglechangeref(ir)) > abs((turnchange-2.0*pi)-anglechangeref(ir))):
				turnchange-= 2*pi
			print("ir = {} turnchange = {} ref= {} ",ir,turnchange*180/pi,anglechangeref[ir]*180/pi)
		turnchangeav+=turnchange
		turna+=turnchange
		print("ir={}, turna = {}, turnchange= ", ir,turna*180/pi,turnchange*180/pi)
		b1+= (ir - incrot)*turna
		b2+=turna
	turnperres= (aa2*b1 - a12*b2)/(a22*a11 - a12*a21)
	turnchangeav=turnchangeav/(nres - 2*incrot -1)
	print("turnchangeav= {} {}",turnchangeav,turnchangeav*180/pi)
	return turnperres



