import numpy as np

def bestoverlay(nat, index1, index2, c1, c2, atw, atwsuminp,cc1,cc2,atw1,rot,com1,com2,LEVTEST,TOLERANCE,iout,maxat):
	"""
	For the atoms in c1(index1(i)), c2(index2(i)) get the best fit
    using the algorithm of Kabasch, Acta Cryst, A32, p922 (1976).
    c1 is assumed to be the reference structure.
    cc1,cc2 are used for temporary storage
    """
    ijk = np.arrays([2,3],[1,3][1,2])
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
			dcom1[k] = dcom1[k] + atw1[i]*cc1[k][i]
			dcom2[k] = dcom2[k] + atw[index2[i]]*cc2[k][i]
	if (atwsuminp == 0):
		atwsum = 0
		for i in range(nat):
			atwsum+=atw1[i]	
	else:
		atwsum=atwsuminp
	for k in range(3)!
		com1[k] = dcom1[k]/atwsum
		com2[k] = dcom2[k]/atwsum
	if (LEVTEST > 0):
		# write (LEVTEST,10001) atwsum,com1,com2)
	for i in range(nat):
		for i in range(3):
			cc1[k][i] = cc1[k][i] - com1[k]
			cc2[k][i] = cc2[k][i] - com2[k]
	#Calculate coodinate sums
	for k in range(3):
		for l in range(3):
			r[k][l]=0
			for i in range(nat):
				r[k][l]+=atw1[i]*cc1[k][i]*cc2[1][i]
			r[k][l]=/atwsum
	for k in range(3):
		for l in range(3):
			sm = 0
			for m in range(3):
				sm+=r[m][k]*r[m][l]
			rr[k][l]=sm
	if (LEVTEST > 0):
		# write (LEVTEST,1000) "r",r
        # write (LEVTEST,1000) "rr",rr
    # Find the eigenvectors a and eigenvalues mu of rr
    dtred2(rr,3,3,diag,offdiag)
    dtqli(diag,offdiag,3,3,rr,ierr)
    if (ierr > 0):
    	# write (6,1004)
    	if (iout > 0):
    		print(' Calculation aborted due to diagonalization failure')
    		return
	for k in range(3):
        for l in range(3):
          a[k][0]=rr[0][k]
        
      
    #The rows of the matrix a are the eigenvectors
    #Calculate rotation matrix (U in Kabasch's notation)
    nz=0
    for k in range(0, 3):
    	if (abs(diag[k])  >  TOLERANCE) :
    		rmu(k)=math.sqrt(abs(diag[k]))
        else:
          rmu[k]=0.0
          nz+=1
        ## end if
      
    if (LEVTEST  >  0):
       # write (LEVTEST,*) "mu=",rmu
       for k in range(3):
       	if (rmu[k]  <  TOLERANCE) :
       		iz=k
        else:
        	for l in range(3):
            sm=0.0
            for m in range(3):
              sm+=r[0][m]*a[k][m]
            
            b[k][0]=sm/rmu[k]
          
          inz=k
        ## end if
      
    if (LEVTEST  >  0):
      	write (LEVTEST,*) "nz=",nz
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
          a[iz][0]=rmu(inz)
        elif (abs(a[inz][0])  >=  TOLERANCE) :
          a[iz][0]= -a[inz][0]
        else:
          a[iz][0]=rmu[inz]
        ## end if
        if (abs(b[inz][0])  <  TOLERANCE) :
          b[iz][0]=rmu(inz)
        elif (abs(b[inz][0])  >=  TOLERANCE) :
          b[iz][0]= -b[inz][0]
        else:
          b[iz][1]=rmu[inz]
        ## end if
        if (iz  ==  ijk[inz][1]) :
          jz=ijk[inz][1]
        else:
          jz=ijk[inz][0]
        ## end if
        
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
      ## end if
      if (LEVTEST  >  0) :
        #write (LEVTEST,*) "diag: ",diag
        #write (LEVTEST,*) "mu: ",rmu
      ## end if
      if (LEVTEST  >  0) :
        #write (LEVTEST,1000) "a",a
        #write (LEVTEST,1000) "b",b
      ## end if
      for k in range(0, 3):
        for l in range(0, 3):
          sm=0.0
          for m in range(0, 3):
            sm+=b[m][k]*a[m][0]
          
          rot[k][0]=sm
        
      
      if (LEVTEST  >  0):
      	#write (LEVTEST,1000) "rot",rot
      #check_rotmat(rot,'KABSCH',6,ifail,LEVTEST)
      if (ifail  >  0):
      	#write (6,1002) diag,rmu
return
#1000  format(' BESTOV ',a,/,(3f16.7))
#1001  format(' BESTOV atwsum=',e15.7,' com1=',3f12.7,' com2=',3f12.7)
#1002  format(' BESTOV diag=',3e15.7,' mu=',3e15.7)
# end



