def trajlimtest(nframe, MAXFRAMES):
	if (nframe > MAXFRAMES):
		raise ValueError(' Calculation stopped since this option is limited to reading {} frames Recompile with larger value of the parameter MAXFRAMES'.format(MAXFRAMES))
		
