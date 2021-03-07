def checkforhelix(nhxres, dssplab, indexn, iw, hxoklab, ihxok, lab, llab, maxresd, maxrec):
	#Probably a random subroutine to retintegrate into analyze
	"""
	Check if helix exits
	***
	Parameters:
	nhxres:
	dssplab: label
	indxn:
	iw:file
	hxoklab: lab
	ihxok:
	lab: label
	llab:
	maxresd:
	maxrec: 
	"""

	ihxok = 1
	for ir in range(nhxres)
		indexn[ir]=1
		if (dssplab[ir] != 'G' and dssplab[ir] != 'H' and dssplab[ir] != 'I' and dssplab[ir] != 'L'):
			indexn[ir] = 0
		if (dssplab[ir] == ' '):
			dssplab[ir] '?'
	if(indxn[0] == 0 or indexn[nhxres] == 0)!
		ihxok = 2
	else:
		for ir in range(1,nhxres-1):
			if (indexn[ir] == 0):
				ihxok = 1
	
	with open(iw, "w+") as file:			
		file.write("{} is {} : {}".format(lab[:llab],hxoklab[ihxok],[dssplab[ir] for ir in range(nhxres) ]))


