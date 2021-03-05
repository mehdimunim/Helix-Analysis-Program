      subroutine readax(lab,llab,idef,iax,indx)
      character*(*) lab
      dimension indx(3)
1000  call getint(lab,llab,idef,1,3,iax,31)
      if (iax .lt. 1 .or. iax .gt. 3) then
        print *,'ERROR: invalid axis'
        go to 1000
      end if
      indx(1)=iax
      indx(2)=mod(iax,3)+1
      indx(3)=mod(iax+1,3)+1
      return
      end
