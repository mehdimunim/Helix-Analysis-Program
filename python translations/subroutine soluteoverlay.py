import numpy as np

def soluteoverlay(isubcrm, ioverlay, nslt, nsegslt, c, cres, ccl, cc2, crmslt0, crmslt, atw, atw1, overlaysds, overlaysh, molsltlim, itemp, idebughx, iw0, maxrsd, maxat):
	# print *,'SOLUTEOVERLAY
	crmshift = np.zeros(3)
	if (isubcrm + ioverlay == 0):
		# No translation or rotation
		cc2 = np.copy(c)
	else:
		cofms(c,crmslt, nslt, atw)
		crmshift = crsmlt - crmslt0
		if (ioverlay == 0):
			#just translation
			for ia in range(nslt):
				cc2[0][ia] = c[0][ia] - crmshift
		elif (ioverlay == 1):
			# overlay of the whole solute
			itemp = [i for i in range(nslt)]
			bestoverlay(nslt, itemp, itemp, cres, c, atw, 0, cc1, cc2, atw1, rot, com1, com2, idebughx, 0.001, iw0, maxat)
			cc2 = shiftmol(c, nslt, com2, -1)
			cc2 = rotate_c(cc2, nslt, rot, 'HELIX', 5)
			cc2 = shiftmol(cc2, nslt, com1, 1)
		overlaysd = sdsumix(nslt, cres, cc2, atw, 0, itemp, devmax, maxat)
		else:
			#Overlay separately each solute molecule
			for iss in range(nsegslt):
				nats = molsltlim[1][iss] - molsltlim[0][iss]+1
				ifat = molsltlim[0][iss]
				itemp = [i for i in range(nats)]
				bestoverlay(nats,itemp,itemp,cres[0][ifat], c[0][ifat], atw[ifat, 0, cc1[0][ifat], cc2[0][ifat], atw1, rot, com1, com2, idebughx, 0.001, iw0, maxat)
				cc2[0][ifat] = shiftmol(c[0][ifat], nats, com2,-1)
				cc2[0][ifat] = rotate_c(cc2[0][ifat], nats, rot , 'HELIXs',6)
				cc2[0][ifat] = shiftmol(cc2[0], ifat, nats, com1, 1)
				overlaysd[iss] = sdsumix(nats, cres[0][ifat], cc2[0][ifat],atw[ifat],0,itemp,devmax,maxat)
				overlaysh[iss] = sqrt(dist2(com1,com2))

	if (ioverlay == 0):
		# write (iw0,1004) crmshift
	elif (ioverlay == 1):
		# write (iw0,1004) crmshift,' RMSD=',overlaysd
	elif (ioverlay == 2):
		# write (iw0,1003) (overlaysh(is),sqrt(overlaysds(is)),is=1,nsegslt)
	elif (isubcrm > 0):
		# write (iw0,1004) crmshift
	return
	# 1003  format(' (Molecular shift, RMSD):',4(' (',f8.2,',',f8.2,')'))
	# 1004  format(' Solute COM shift=',3f8.2,a,f8.2)
