      subroutine check_rotmat(rot,lab,llab,ifail,iverb)
      dimension rot(3,3)
      character*(*) lab
      common /rotwarn/ nrotwarn,nrottest,devmax
      dimension u(3,3)
      data u /1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0/
c     Check the rotation matrix
      ifail=0
      nrottest=nrottest+1
      do i=1,3
        do j=1,3
          sum=0.0
          do k=1,3
            sum=sum+rot(i,k)*rot(j,k)
          end do
          if (abs(sum-u(i,j)) .gt. 0.99) then
            ifail=1
            write (6,1000) 'PROGRAM ERROR',lab(1:llab),i,j,sum
c           print *,'STOPPING'
c           STOP
          else if (abs(sum-u(i,j)) .gt. 0.1) then
            if (nrotwarn .lt. 00)
     -        write (6,1000) 'WARNING',lab(1:llab),i,j,sum
            nrotwarn=nrotwarn+1
            if (devmax .lt. abs(sum-u(i,j))) devmax=abs(sum-u(i,j))
          end if
        end do
      end do
      if (ifail .gt. 0 .or. iverb .gt. 0) write (6,1001) rot
      return
1000  format(1x,a,': rotation matrix error in call ',a,/,
     -  ' SUM rot(',i1,',k)*rot(',i1,',k)=',e12.5)
1001  format(' The rotation matrix in rotate_c:',/,(3e13.5))
      end
