def compared2(es,ex,d2min,hl,hlax,a1,a2,aes):
	"""
	
	Compare two vectors
	NOTE: in comment in simulaid.f

	***
	Parameters:


	"""
	if (aes > -h1 and aes <= h1):
		d2 = dist2(es,ex)
		if (d2min > d2):
			d2min = d2
			a1 = hlax
			a2 = aes
	return d2