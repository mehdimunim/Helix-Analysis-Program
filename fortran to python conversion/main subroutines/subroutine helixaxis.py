import numpy as np
import math


def helixaxis(nframe, max0, c,nslt,iw0,calph,axisdir,axisini,axisend,helixcent,perpvec,camod,anglechangeref, rn,axfact,axtol,rot,rms,helixlen,angles,decidebend,nup,ndown,nrun,nnear,rcirc,turnperres,anglesn,irot,incrot,nrep,irefang,nreshx,icaahx,ihx,nhxres,idebughx,radtodeg,maxhlx):
    """
    This method return a condensed result summarizing the info about the helix
    
    ***
    Parameters:
    nframe:
    max0:
    c: input structure
    nslt: number of solvent molecules
    iw0:
    calph: carbone alpha
    perpvec:
    camod:
    anglechangeref:
    rn: 
    axfact: 
    axtol:
    rot: 
    decidebend:
    irot:
    incrot:
    nrep:
    irefang:
    nreshx:
    icaahx:
    ihx:
    nhxres:
    idebughx:
    radtodeg:
    maxhlx:
    
    Returns:
    res: a condensed result

    """
    MAXFRAMES=50000
    MAXCOPY=600
    res = np.zeros(2,MAXFRAMES,MAXCOPY)
    decide = ['Too short','Bent','Random','Alternate','Too long ','']
    print("HELIXAXIS maxhlx,nreshx,nrep = ",maxhlx,nreshx,nrep)
    iprintkahn=0
    it=0
    for ir in range(nreshx):
    	for k in range(3):
    		calph[k][ir]=c[k][icaahx[ir]]
    if (irot  ==  1) :
    	for ir in range(nreshx):
    		calph[1][ir] = dsmatvec(rot,calph[0][ir])
    message='Helix'
    rms, axisdir, axisini, axisend = kahn(calph,nreshx, True ,axisdir,axisini,axisend,rms,iprintkahn, message,maxhlx)
    rmsa=rms
    if (idebughx  >  0) :
        print("Helix axis direction =",axisdir)
        print("Helix axis position = ", axisini)
        print("Helix end position = ", axisend)
        print("Helix quality rms = ",rms)
        print("# of HX res = ",nreshx)
        for ir in range(nreshx):
            print(icaahx[ir])
    if (nrep  <=  1) :
    	for ir in range(nreshx):
          perpvec[0][ir], camod[0][ir] = calcperp(axisini,axisdir,calph[0][ir],it)
          ibtyp, nup,ndown,nrun,nnear, rcirc,helixcent,rn = checkbend(calph,axisdir,camod,axfact,perpvec, nreshx,axtol,idebughx)
    if (nrep  <=  1) :
    	for k in range(3):
            try:
    		  anglesn[k]=dacoscheck(rn[k])
            except:
                print("ERROR in dacoscheck angle=",rn[k])
        helixlen=math.sqrt(ddistsq(axisini,axisend))
        for k in range(3):
            try:
              angles[k]=dacoscheck(axisdir[k])
            except:
                print("ERROR in dacoscheck angle=",axisdir[k])
        turnperres = calcturnperres(nreshx,incrot,perpvec,axisdir,anglechangeref,irefang,maxhlx)
        nframes=max0[0][nframe]
        incr=(ihx-1)*nhxres
        
        try:
            trajlimtest(nframe,MAXFRAMES)
        except:
            print("ERROR with number of frames")

        res[0][nframes][incr+7]=helixlen
        res[1][nframes][incr+7]=rcirc
        for k in range(3):
          res[0][nframes][incr+k]=math.cos(angles[k])
          res[1][nframes][incr+k]=math.sin(angles[k])
        res[0][nframes][incr+6]=math.cos(turnperres)
        res[1][nframes][incr+6]=math.sin(turnperres)
        res[0][nframes][incr+13]=axisdir[0]
        res[1][nframes][incr+13]=axisdir[1]
        res[0][nframes][incr+14]=axisdir[2]
        res[1][nframes][incr+14]=helixcent[0]
        res[0][nframes][incr+15]=helixcent[1]
        res[1][nframes][incr+15]=helixcent[2]

        decidebend=decide[ibtyp]

        if (iw0  >  0):
        	printhelix(axisini,axisend,helixcent,rms,helixlen,axisdir,angles,decidebend,nup,ndown,nrun,nnear,rcirc,turnperres,anglesn,ihx,radtodeg)
    return res