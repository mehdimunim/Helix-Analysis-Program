import numpy as np

def cofms(csa, aw):
	"""
	Compute center of mass from global atomic coordinates
			SUM(csax[i]*aw[i]) + SUM(csay[i]*aw[i]) + SUM(csaz[i]*aw[i])
	rmass = -----------------------------------------------------------
			                SUM(aw[i])
	***
	Parameters:
	csa: input coordinates
	aw: atom weights
	"""
	#number of atoms
	natm = np.size(csa[0])
	#calculate rmass numerator
	for i in range(3):
		rmass[i]=sum([aw[j]*csa[i][j] for j in range(natm)])
	#calculate rmass denominator, wmol the total weight of molecule
	wmol = sum(aw)
	#if wmol is not defined, make it not account
	if (wmol == 0):
		wmol=1
	return rmass/wmol	