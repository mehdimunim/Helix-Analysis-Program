      subroutine genrot(rot,pi,iax,angle)
      dimension rot(3,3)
      character*1 ans
      character*6 rq
      print *,'Rotating the structure (if required)'
      call quiz(ans,iax,' ',' ',0,'rotation type',13,0,5,6,0)
c     print *,'ans=',ans
      if (ans .eq. 'i') then
        iax=0
      else if (ans .eq. 'q') then
        iax=-1
        return
      end if
      if (iax .gt. 0) then
        call getreal(
     -    'Angle (+: clockwise, -: counterclockwise, viewed from +'
     -    //ans//' axis)',62,999999.0,angle,0,0)
        angler=angle*(pi/180.0)
        iz=iax
        ix=mod(iz,3)+1
        iy=mod(iz+1,3)+1
        call unitmat(rot)
c       print *,'ix,iy,iz=',ix,iy,iz
        rot(ix,ix)=cos(angler)
        rot(ix,iy)=sin(angler)
        rot(iy,ix)=-rot(ix,iy)
        rot(iy,iy)=rot(ix,ix)
      else
c       Input rotation matrix
        do i=1,3
          do j=1,3
            write (rq,1003) i,j
            call getreal(rq,6,999999.0,rot(i,j),0,00)
          end do
        end do
      end if
c      do i=1,3
c        write (6,7777) i,(rot(j,i),j=1,3)
c7777    format(i3,' rot i =',3f7.4)
c      end do
      return
1003  format('R(',i1,',',i2,')')
      end
