import numpy as np

def prokinkcalcla(nrepc, c, nslt, bend, wobble, facshift, rmsb, rmsa, iflatproline, ix5, radtodeg, maxhlx):
	calph = np.zeros(3,maxhlx)
	axisdirb = np.zeros(3)
	axisinib = np.zeros(3)
	axisendb = np.zeros(3)
	aixsdira = np.zeros(3)
	axisinia = np.zeros(3)
	axiseenda = np.zeros(3)
	rms = 0
	MAXHX = 50
	#ProKink variable
	ddot = 0
	dmag = 0
	cosa = 0
	pro_alphaC = 0
	orig = 0
	origa = 0
	orig3 = 0
	orig4 = 0
	perpvec = 0
	perpveca = 0
	perpvec3 = 0
	perpvec4 = 0
	perpvec34 = 0
	enda = 0
	bmin3C = 0 
    bmin4C = 0
    zax = 0
    xx = 0
    prolinering = 0
    p0 = 0
    ringnorm = 0
    axfact = 0
    rr = 0
    pro_alphaC = np.zeros(3)
    orig = np.zeros(3)
    origa = np.zeros(3)
    orig4 = np.zeros(3)
    perpvec = np.zeros(3)
    perveca = np.zeros(3)
    perpvec3 = np.zeros(3)
    perpvec4 = np.zeros(3)
    perpvec34 = np.zeros(3)
    bmin3C = np.zeros(3)
    bmin4C = np.zeros(3)
    enda = np.zeros(3)
    zax = np.zeros(3)
    xx = np.zeros(3)
    prolinering = np.zeros(3,5)
    p0 = np.zeros(3)
    ringnorm = np.zeros(3)
    axfact = np.zeros(5)
    iprintkahn = 0
    for k in range(3):
    	calph[k][0]=c[k][icapr]
    for ir in range(nra):
    	for k in range(3):
    		calph[k][ir +1] = c[k][icaa[ir+1]]
    message ="Helix after the kink"
    kahn(calph,nra+1,.true.,axisdira,axisinia,axisenda,rms,iprintkahn, message,MAXHX)
    rmsa = rms
    if (iprintpk > 1):
    	print("Helix after: axis",axisdira)
    	print("Helix after: init",axisinia)
    	print("Helix after: end",axisenda)
    	print("RMS =",rms)
    for ir in range(nrb):
    	for k in range(3):
    		# calph[k][ir+1]=c[k][icab[ir+1]]
    		calph[k][ir]=c[k][icab[ir]]
    message = "Helix before the kink"
    kahn(calph,nrb,.true.,axisdirb,axisinib,axisendb,rms, iprintkahn, message,MAXHX)
    rmsb = rms
    if (iprintpk > 1):
    	print("Helix before: axis=",axisdirb)
    	print("Helix before: init=",axisinib)
    	print("Helix before: end=",axisendb)
    	print("RMS= ",rms)
    iac = icapr
    im3 = icab[2]
    im4 = icab[3]
    for k in range(3):
    	pro_alphaC[k]=c[k][icapr]
    	bmin3C[k]=c[k][im3]
    	bmin4C[k]=c[k][im4]
    if (iflatproline == 1):
    for i in range(5):
    	for k in range(3):
    		prolinering[k][i] = c[k][ix5[i]]
    fitpoints(prolinering,5,3,p0,ringnorm,axfact,iprintpk)
    xx = pro_alphaC - p0
    rr = np.dot(xx,ringnorm)
    for k in range(3):
    	pro_alphaC[k] +=rr*ringnorm[k]
    # New code for proline kink calculation
    # Bend angle: from the scalar product of the two axis vectors
    # Both axis vectors point away from the proline
    coa = - np.dot(axisdira,axisdirb)
    bend = dacoscheck(cosa,ccc,1,6,"PROK-B")*radtodeg
   	# Wobble angle: from the scalar product of the two normals to the 
   	# before helix axis (from the C-alpha of Proline and the end of the after helix) 
   	calcperp(axisendb,axisdirb,pro_alphC,orig,perpvec,iprintpk)
   	enda = orig + axisdira
   	calcperp(axisinib,axisdirb,enda,origa,perpveca,iprintpk)
   	xx = origa - enda
   	cc = np.linalg.norm(xx)
   	if (cc < 0.01):
   		print ("Wobble is suspect")
   	cosa = np.dot(perpvec,perpveca)
   	wobble = dacoscheck(cosa,ccc,1,6,'PROK-W')*radtodeg
   	#Establish sign
   	dcross(axisdirb,perpvec,zax)
   	if (np.dot(zax,perpveca)>0):
   		wobble = -wobble
   	calcperp(axisinib,axisdirb,bmin3C,orig3,perpvec3,iprintpk)
   	calcperp(axisinib,axisdirb,bmin4C,orig4,perpvec4,iprintpk)
   	perpvec = dvnorm(perpvec)
   	perpvec3 = dvnorm(perpvec3)
   	perpvec4 = dvnorm(perpvec4)
   	for k in range(3):
   		perpvec34[k] = (perpvec3[k] + perpvec4[k])/2
   	perpvec34 = dvnorm(perpvec34)
   	cosa = np.dot(perpvec,perpvec34)
   	facshift = dacoscheck(cosa,ccc,1,6,"PROK-FS")*radtodeg
   	# Establish sign
   	cosa3 = np.dot(perpvec,perpvec3)
   	cosa4 = np.dot(perpvec,perpvec4)
   	if (cosa3 > cosa4):
   		faceshift = -faceshift
   	if (nrepc < 0):
   		return none
   	return faceshift






