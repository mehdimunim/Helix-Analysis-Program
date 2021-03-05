def printhelix(iw0, axisini, axisend, cent, rms, helixlen, axisdir, angles, decidebend, nup, ndown, nrun, nnrear, rcirc, turnperres, anglesn, ihx, radtodeg):
	"""
	Print an helix
	***
	Parameters:
	iw0: file
	axisini: start of the axis
	axisend: end of the axis
	cent:
	rms: root mean square
	helixlen: length of the helix
	axisdir: axis vector
	angles: angle between X,Y, Z and the axis of the helix
	decidebend: Shape of the helix
	nup:
	ndown:
	nrun:
	nnear:
	rcirc:
	turnperres:TPR/ Turn Per Residue
	anglesn:
	ihx: Number of the helix
	radtodeg: Conversion factor from radian to degrees
	"""
	with open(iw0,"w+") as file:
		file.write("HX# {} S= {} E= {} Len= {} D= {}".format(ihx,axisini,axisend,helixlen,axisdir))
		file.write("D-X,D-Y,D-Z angles {} RMS= {} Shape = {}".format([radtodeg*angles[k] for k in range(3)],rms,decidebend))
		file.write("Nup/dn= {} Ncross = {} Nax = {} Rc= {} TPR = {} C= {} N-X, N-Y, N-Z angles ={}".format(nup,ndown,nrun-1,nnear,rcirc,radtodeg*turnperres,cent,[radtodeg*anglesn[k] for k in range(3)]))

