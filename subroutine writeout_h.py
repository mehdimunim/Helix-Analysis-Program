def writeout_h(dir, ip, fp, rms, message, iprint):
	"""
	Write header???
	***
	Parameters:
	dir: direction vector
	ip: initial point
	fp: final point
	rms: root mean square
	"""
	kprint = mod(iprint,10)
	# don't waste time
	if (kprint == 0):
		return
	if (message[:5] != ""):
		print(message)
	print("D ",dir)
	print("I",ip)
	print("F,",fp)
	print("RMS deviation:",rms)