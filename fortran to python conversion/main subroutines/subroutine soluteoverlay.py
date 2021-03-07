import numpy as np

def soluteoverlay(isubcrm, ioverlay, nslt, nsegslt, c, cres, crmslt0, crmslt, atw, atw1, overlaysds, overlaysh, molsltlim, itemp, idebughx, iw0, maxrsd, maxat):
	"""
	
	Translate and rotate the solute

	***
	Parameters:
	isubcrm:
	ioverlay:
	nslt: number of solvant molecules
	nsegslt:
	maxat:
	crmslt0:

	Returns:

	"""
	crmshift = np.zeros(3)
	if (isubcrm + ioverlay == 0):
		# No translation or rotation
		cc2 = np.copy(c)
	else:
		crmslt = cofms(c, nslt, atw)
		#shift the center of mass
		crmshift = crsmlt - crmslt0
		if (ioverlay == 0):
			#just translation
			for ia in range(nslt):
				cc2[0][ia] = c[0][ia] - crmshift
		elif (ioverlay == 1):
			# overlay of the whole solute
			itemp = [i for i in range(nslt)]
			com1, com2, rot = bestoverlay(nslt, itemp, itemp, cres, c, atw, 0, cc1, cc2, idebughx, 0.001, iw0, maxat)
			cc2 = shiftmol(c, nslt, com2, -1)
			cc2 = rotate_c(cc2, nslt, rot)
			cc2 = shiftmol(cc2, nslt, com1, 1)
			overlaysd = sdsumix(nslt, cres, cc2, atw, 0, itemp, devmax, maxat)
		else:
			#Overlay separately each solute molecule
			for iss in range(nsegslt):
				nats = molsltlim[1][iss] - molsltlim[0][iss]+1
				ifat = molsltlim[0][iss]
				itemp = [i for i in range(nats)]
				com1, com2, rot = bestoverlay(nats,itemp,itemp,cres[0][ifat], c[0][ifat], atw[ifat], 0, cc1[0][ifat], cc2[0][ifat], atw1, rot, com1, com2, idebughx, 0.001, iw0, maxat)
				cc2[0][ifat] = shiftmol(c[0][ifat], nats, com2,-1)
				cc2[0][ifat] = rotate_c(cc2[0][ifat], nats, rot , 'HELIXs',6)
				cc2[0][ifat] = shiftmol(cc2[0], ifat, nats, com1, 1)
				overlaysd[iss] = sdsumix(nats, cres[0][ifat], cc2[0][ifat],atw[ifat],0,itemp,devmax,maxat)
				overlaysh[iss] = math.sqrt(dist2(com1,com2))

	if (ioverlay == 0):
		print("Solute COM shift =",crmshift)
	elif (ioverlay == 1):
		print("Solute COM shift= {} RMSD = {}".format(crmshift,overlaysd))
	elif (ioverlay == 2):
		for iss in range(nsegslt):
			print("(Molecular shift, RMSD: ",overlaysh[iss],math.sqrt(overlaysds[iss]))
	elif (isubcrm > 0):
		print(" Solute COM shift=",crmshift)
	return crmshift, overlaysd
