      subroutine transp(rot,N)
c     /*- transpose an NxN matrix rot -*/
      real rot(N,N),tmp
      do i=1,N
         do j=i+1,N
            tmp=rot(j,i)
            rot(j,i)=rot(i,j)
            rot(i,j)=tmp
         end do
      end do
      return
      end
