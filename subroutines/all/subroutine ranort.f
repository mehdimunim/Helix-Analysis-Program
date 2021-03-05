      subroutine ranort(r,i2dopt)
c#    MMC routine 097 lstmod: 10/10/90
c*****Generate an orientation matrix from the uniform distribution
c     from the 3 Euler angles (i2dopt=0) or a random orientation, corresponding
c     to a random rotation around axis i2dopt (i2dopt .gt. 0)
      dimension r(3,3),rn(3)
      if (i2dopt .eq. 0) then
c       eq(4-47) of Goldstein (r=a)
        call randpx(3,rn)
        fir=rn(1)*3.141592*2.0
        Psr=rn(3)*3.141592*2.0
        Cth=-1.0+2.0*Rn(2)
        sth=sqrt(1.0-Cth*cth)
        sfi=sin(fir)
        cfi=cos(fir)
        sps=sin(psr)
        cps=cos(psr)
        r(1,1)=cps*cfi-cth*sfi*sps
        r(2,1)=-sps*cfi-cth*sfi*cps
        r(3,1)=sth*sfi
        r(1,2)=cps*sfi+cth*cfi*sps
        r(2,2)=-sps*sfi+cth*cfi*cps
        r(3,2)=-sth*cfi
        r(1,3)=sth*sps
        r(2,3)=sth*cps
        r(3,3)=cth
      else
        call randpx(1,rn)
        fir=rn(1)*3.141592*2.0
        sfi=sin(fir)
        cfi=cos(fir)
        call zeroit(r,9)
        r(i2dopt,i2dopt)=1.0
        ix=mod(i2dopt,3)+1
        iy=mod(i2dopt+1,3)+1
        r(ix,ix)=cfi
        r(ix,iy)=sfi
        r(iy,ix)=-sfi
        r(iy,iy)=cfi
      end if
      return
      end
