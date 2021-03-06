      subroutine euler(r,fi,ps,th)
c#    MMC routine 097 lstmod: 10/10/90
c*****Generate an orientation matrix from the 3 euler angles
      dimension r(3,3)
c     this subroutine prepares a rotation matrix
c     eq(4-47) of Goldstein (r=a)
      sfi=sin(fi)
      cfi=cos(fi)
      sps=sin(ps)
      cps=cos(ps)
      sth=sin(th)
      cth=cos(th)
      r(1,1)=cps*cfi-cth*sfi*sps
      r(2,1)=-sps*cfi-cth*sfi*cps
      r(3,1)=sth*sfi
      r(1,2)=cps*sfi+cth*cfi*sps
      r(2,2)=-sps*sfi+cth*cfi*cps
      r(3,2)=-sth*cfi
      r(1,3)=sth*sps
      r(2,3)=sth*cps
      r(3,3)=cth
      return
      end
