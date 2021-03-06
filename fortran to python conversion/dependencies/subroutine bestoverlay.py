import numpy as np

def bestoverlay(nat, index1, index2, c1, c2, atw, atwsuminp, atw1, LEVTEST,TOLERANCE,maxat):
	"""
  Compare c1 and c2 by calculating the rotation matrix that
  minimizes the RMSD (root mean square deviation)
	For the atoms in c1(index1(i)), c2(index2(i)) get the best fit
  using the algorithm of Kabsch, Acta Cryst, A32, p922 (1976).
  c1 is assumed to be the reference structure.
  cc1,cc2 are used for temporary storage
  ***
  Parameters:
  nat: number of atom
  index1:
  index2:
  c1: reference structure
  c2: another structure
  atw: atomic weight
  atwsuminp:
  atw1: atom weight 1
  LEVTEST:
  TOLERANCE:
  maxat: 

  Returns:
  rot: optimized rotation matrix
  """
  com1 = np.zeros(3)
  com2 = np.zeros(3)
  cc1 = np.zeros(3,nat)
  cc2 = np.zeros(3,nat)
  print("BESTOVERLAY LEVTEST,maxat,nat=",LEVTEST,maxat,nat)
  ijk = np.arrays([2,3],[1,3][1,2])
  devmax = 0
	if (nat == 1):
		return
	# Move both sets into COM frame
	for k in range(3):
		dcom1[k]= 0
		dcom2[k]=0
	for i in range(nat):
		atw1[i] = atw[index1[i]]
		for k in range(3):
			cc1[k][i] = c1[k][index1[i]]
			cc2[k][i] = c2[k][index2[i]]
			dcom1[k] += atw1[i]*cc1[k][i]
			dcom2[k] += atw[index2[i]]*cc2[k][i]
	if (atwsuminp == 0):
		atwsum = 0
		for i in range(nat):
			atwsum += atw1[i]	
	else:
		atwsum=atwsuminp
	for k in range(3)!
		com1[k] = dcom1[k]/atwsum
		com2[k] = dcom2[k]/atwsum
	if (LEVTEST > 0):
    print("BESTOV atwsum={} com1={} com2={}".format(atwsum,com1,com2))
	for i in range(nat):
		for i in range(3):
			cc1[k][i] = cc1[k][i] - com1[k]
			cc2[k][i] = cc2[k][i] - com2[k]
	#Calculate coodinate sums
	for k in range(3):
		for l in range(3):
			r[k][l]=0
			for i in range(nat):
				r[k][l]+=atw1[i]*cc1[k][i]*cc2[l][i]
			r[k][l]/=atwsum
	for k in range(3):
		for l in range(3):
			sm = 0
			for m in range(3):
				sm+=r[m][k]*r[m][l]
			rr[k][l]=sm
	if (LEVTEST > 0):
    print("BESTOV r",r)
    print("BESTOV rr",r)
	# Find the eigenvectors a and eigenvalues mu of rr
  diag,offdiag,rr = dtred2(rr,3,3,diag,offdiag)
  diag,offdiag,rr, ierr = dtqli(diag,offdiag,3,3,rr,ierr)
  # to be modified as in checkbend
  if (ierr > 0):
    # write (6,1004)
    if (iout > 0):
    	print(' Calculation aborted due to diagonalization failure')
	for k in range(3):
        for l in range(3):
          a[k][l]=rr[l][k]
  #The rows of the matrix a are the eigenvectors
  #Calculate rotation matrix (U in Kabsch's notation)
  nz=0
  for k in range(0, 3):
    if (abs(diag[k])  >  TOLERANCE) :
    	rmu[k]=math.sqrt(abs(diag[k]))
    else:
      rmu[k]=0.0
      nz+=1
  if (LEVTEST  >  0):
    print("mu=",rmu)
  for k in range(3):
    if (rmu[k]  <  TOLERANCE) :
      iz=k
    else:
      for l in range(3):
        sm=0.0
        for m in range(3):
          sm+=r[l][m]*a[k][m]
        b[k][l]=sm/rmu[k]
      inz=k
  if (LEVTEST  >  0):
    print("nz = ",nz)
  if (nz  ==  1) :
    # Planar molecule: a(iz)=a(inz)xa(jnz), same for b
    a[iz][0]=a[ijk[iz][0],1]*a[ijk[iz][1],2]- a[ijk[iz][0],2]*a[ijk[iz][1],1]

    a[iz][1]=a[ijk[iz][0],2]*a[ijk[iz][1],0]- a[ijk[iz][0],0]*a[ijk[iz][1],2]

    a[iz][2]=a[ijk[iz][0],0]*a[ijk[iz][1],1]- a[ijk[iz][0],1]*a[ijk[iz][1],0]

    b[iz][0]=b[ijk[iz][0],1]*b[ijk[iz][1],2]- b[ijk[iz][0],2]*b[ijk[iz][1],1]

    b[iz][1]=b[ijk[iz][0],2]*b[ijk[iz][1],0]- b[ijk[iz][0],0]*b[ijk[iz][1],2]

    b[iz][2]=b[ijk[iz][0],0]*b[ijk[iz][1],1]- b[ijk[iz][0],1]*b[ijk[iz][1],0]
  elif (nz  ==  2) :
    # Linear molecules
    # First set the iz-th components to nonparalel to the inz-th
    for k in range(3):
      a[iz][k]=0.0
      b[iz][k]=0.0
    if (abs(a[inz][0])  <  TOLERANCE) :
      a[iz][0]=rmu[inz]
    elif (abs(a[inz][0])  >=  TOLERANCE) :
      a[iz][0]= -a[inz][0]
    else:
      a[iz][1]=rmu[inz]
    if (abs(b[inz][0])  <  TOLERANCE) :
      b[iz][0]=rmu(inz)
    elif (abs(b[inz][0])  >=  TOLERANCE) :
      b[iz][0]= -b[inz][0]
    else:
      b[iz][1]=rmu[inz]
    if (iz  ==  ijk[inz][0]) :
      jz=ijk[inz][1]
    else:
      jz=ijk[inz][0]
    # a(jz)=a(inz)xa(iz)
        
    a[jz][0]=a[ijk[jz][0],1]*a[ijk[jz][1],2]- a[ijk[jz][0],2]*a[ijk[jz][1],1]

    a[jz][1]=a[ijk[jz][0],2]*a[ijk[jz][1],0]- a[ijk[jz][0],0]*a[ijk[jz][1],2]

    a[jz][2]=a[ijk[jz][0],0]*a[ijk[jz][1],1]- a[ijk[jz][0],1]*a[ijk[jz][1],0]

    b[jz][0]=b[ijk[jz][0],1]*b[ijk[jz][1],2]- b[ijk[jz][0],2]*b[ijk[jz][1],1]

    b[jz][1]=b[ijk[jz][0],2]*b[ijk[jz][1],0]- b[ijk[iz][0],0]*b[ijk[jz][1],2]

    b[jz][2]=b[ijk[jz][0],0]*b[ijk[jz][1],1]- b[ijk[jz][0],1]*b[ijk[jz][1],0]
        
    #a(iz)=a(inz)xa(jz)
    a[iz][0]=a[ijk[iz][0],1]*a[ijk[iz][1],2]- a[ijk[iz][0],2]*a[ijk[iz][1],1]

    a[iz][1]=a[ijk[iz][0],2]*a[ijk[iz][1],0]- a[ijk[iz][0],0]*a[ijk[iz][1],2]

    a[iz][2]=a[ijk[iz][0],0]*a[ijk[iz][1],1]- a[ijk[iz][0],1]*a[ijk[iz][1],0]

    b[iz][0]=b[ijk[iz][0],1]*b[ijk[iz][1],2]- b[ijk[iz][0],2]*b[ijk[iz][1],1]

    b[iz][1]=b[ijk[iz][0],2]*b[ijk[iz][1],0]- b[ijk[iz][0],0]*b[ijk[iz][1],2]

    b[iz][2]=b[ijk[iz][0],0]*b[ijk[iz][1],1]- b[ijk[iz][0],1]*b[ijk[iz][1],0]
    
  if (LEVTEST  >  0) :
    print("diag:",diag)
    print("mu: ",rmu)
  if (LEVTEST  >  0) :
    print("BESTOV a:",a)
    print("BESTOV b: ",b)
  for k in range(3):
    for l in range(3):
      sm=0.0
      for m in range(0, 3):
        sm+=b[m][k]*a[m][l]
      rot[k][0]=sm
  if (LEVTEST  >  0):
    print("BESTOV rot",rot)
  try:
    check_rotmat(rot)
  except:
    print("BESTOV failed rot=",rot)
  print("BESTOV diag={} mu={}".format(diag,rmu))
  return rot