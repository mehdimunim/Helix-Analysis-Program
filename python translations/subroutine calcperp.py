import numpy as np

def calcperp(start,dir,fro,orig, perpvec,itest):
	"""
	For a line from start in the direction dir, calculate the normal to it
	from the point fro. The normal meets the line at orig and its direction
	is perpvec

	"""
	dsx = np.zeros(3)
	xx = np.zeros(3)
	dsx = fro - start
	dsx = np.dot(dir,dsx)
	for k in range(3):
		perpvec[k] = dir[k]*ddsx - dsx[k]
	perpvec = ddnorm(perpvec)
	perpvecfac = np.dot(perpvec,dsx)
	for k in range(3):
		orig[k] = fro[k] - perpvecfac*perpvec[k]
	if (itest > 0):
		print('Start',start)
        print('Dir {} Magn = {}',dir,np.linalg.norm(dir))
        print('From',fro)
        print('Orig',orig)
        print('Perpvec {} Magn = {}',perpvec,np.linalg.norm(dir))
        if (itest > 1):
        	cc = np.dot(dir,perpvec)
        	print("Dir . Perpevec (->0)=",cc)
        	xx = start - orig
        	cc = np.dot(dir,xx)/np.linalgn.norm(xx)
        	print("Dir . (orig-start)",cc)
   	return cc

