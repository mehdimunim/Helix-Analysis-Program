def helixaxis(c,nslt,iw0,calph,axisdir,axisini,axisend,helixcent,perpvec,camod,anglechangeref,circ,rn,axfact,axtol,rot,rms,helixlen,angles,decidebend,nup,ndown,nrun,nnear,rcirc,turnperres,anglesn,irot,incrot,nrep,irefang,nreshx,icaahx,ihx,nhxres,idebughx,radtodeg,pi,maxhlx):
    MAXFRAMES=50000
    MAXCOPY=600
    #common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,increment,increment2,res(2,MAXFRAMES,MAXCOPY),x0(MAXCOPY),y0(MAXCOPY),nxselres,ixselres(MAXCOPY)
    #real*8 ddistsq
    #character*9 decide(6)
    decide = ['Too short','Bent','Random','Alternate','Too long ','']
    #print *,'HELIXAXIS maxhlx,nreshx,nrep=',maxhlx,nreshx,nrep
    iprintkahn=0
    it=0
    for ir in range(0, nreshx):
    	for k in range(0, 3):
    		calph[k][ir]=c[k][icaahx[ir]]
    if (irot  ==  1) :
    	for ir in range(0, nreshx):
    		dsmatvec(rot,calph[0][ir],calph[1][ir])
    message='Helix'
    kahn(calph,nreshx, True ,axisdir,axisini,axisend,rms,iprintkahn, message,maxhlx)
    rmsa=rms
    if (idebughx  >  0) :
    	# write (77,7011) 'axis direction',axisdir
        # write (77,7011) 'initial position',axisini
        # write (77,7011) '# end position',axis# end
        # write (77,*) 'Helix quality rms=',rms
        # write (77,7012) nreshx,(icaahx(ir),ir=1,nreshx)
    if (nrep  <=  1) :
    	for ir in range(0, nreshx):
          calcperp(axisini,axisdir,calph[1][ir],camod[1][ir], perpvec[1][ir],it)
          checkbend(calph,axisdir,camod,axfact,perpvec,helixcent,nreshx,nup,ndown,nrun,nnear,axtol,rcirc,circ,rn,ibtyp,idebughx)
    if (nrep  <=  1) :
    	for k in range(0, 3):
    		anglesn[k]=dacoscheck(rn[k],ccc,1,6,'HELIXNORMAL')
        helixlen=math.sqrt(ddistsq(axisini,axis))
        for k in range(0, 3):
          angles[k]=dacoscheck(axisdir(k),ccc,1,6,'HELIXAXIS')
        calcturnperres(turnperres,nreshx,incrot,perpvec,axisdir,anglechangeref,irefang,pi,maxhlx)
        nframes=max0[1][nframe]
        incr=(ihx-1)*nhxres
        trajlimtest(nframe,MAXFRAMES)
        res(1,nframes,incr+7)=helixlen
        res(2,nframes,incr+7)=rcirc
        for k in range(0, 3):
          res(1,nframes,incr+k)=cos(angles(k))
          res(2,nframes,incr+k)=sin(angles(k))
        res(1,nframes,incr+6)=cos(turnperres)
        res(2,nframes,incr+6)=sin(turnperres)
        res(1,nframes,incr+13)=axisdir(1)
        res(2,nframes,incr+13)=axisdir(2)
        res(1,nframes,incr+14)=axisdir(3)
        res(2,nframes,incr+14)=helixcent(1)
        res(1,nframes,incr+15)=helixcent(2)
        res(2,nframes,incr+15)=helixcent(3)
        decidebend=decide[ibtyp]
        if (iw0  >  0):
        	printhelix(iw0,axisini,axisend,helixcent,rms,helixlen,axisdir,angles,decidebend,nup,ndown,nrun,nnear,rcirc,turnperres,anglesn,ihx,radtodeg)
    return
	# 7011  format(' Helix ',a,'=',3f10.4)
	# 7012  format(' # of HX res=',i2,(' icaa:',20i5))
