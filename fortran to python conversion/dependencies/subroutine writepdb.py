def writepdb(iout,c,ia,ir,atnam,resnam,segnam,frocc,bfac):
	"""
	Write file iout in PDB format
	***
	Parameters
	iout: output file
	ia: 
	atnam: atom name
	resnam:
	segnam: segment id
	frocc:
	bfac:
	"""
	with open(iout,"w+") as f:
		f.write("ATOM {} {} {} {} {} {} {} {}".format(ia,atnam,resnam,segnam,ir,c,frocc,bfac))
