      subroutine fixup(rc,co,cn,n,aw)
      dimension rc(3,n),co(3,n),cn(3,n),aw(n)
c     The molecule co will be replaced by a molecule having the same
c     internal geometry as rc with the com's coinciding and the best
c     orientational overlap
      dimension como(3),a(3,3),b(3,3),orient(3,3)
c     Obtain  orientation of co
c     print *,'fixup n=',n
c     write (6,1002) co
c1002  format(' co=',/,(3f10.5))
      call cofms(co,como,n,aw)
      call trnsfr(a,rc,9)
      do i=1,3
        do j=1,3
          b(i,j)=co(i,j)-como(i)
        end do
      end do
      call ormat(orient,a,b,n,1,6)
      call rotate_c(rc,n,orient,cn,'FIXUP',5)
      call shiftmol(cn,n,como,cn,1.0)
      return
      end
