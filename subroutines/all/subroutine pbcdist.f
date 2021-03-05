      subroutine pbcdist(c1,c2,i1,i2,cell,ncell,iw,nframe,img,r,dist)
      dimension c1(3),c2(3),cell(3,ncell),r(3)
      call arrdiff(c2,c1,r,3)
      call genimdist(r,cell,1,ncell,img,d2)
      dist=sqrt(d2)
      if (img .gt. 1) then
        do k=1,3
          r(k)=r(k)-cell(k,img)
        end do
        if (iw .lt. -1) write (-iw,1002) nframe,img
      end if
      if (iw .eq. 0) write (6,1000) i1,i2,dist,r,img
      if (iw .gt. 0) write (iw,1001) nframe,i1,i2,dist,r,img
      return
1000  format(' D(',i5,'-',i5,')=',f9.3,' r=',3f9.4,
     -  ' A, PBC image=',i2)
1001  format(' N=',i6,' D(',i5,'-',i5,')=',f9.3,' r=',3f9.4,
     -  ' A, PBC image=',i2)
1002  format(' PBC reset at frame # ',i6,' PBC image=',i2)
      end
