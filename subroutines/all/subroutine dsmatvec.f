      subroutine dsmatvec(rot,din,dout)
      real*8 din(3),dout(3),d(3)
      dimension rot(3,3)
      do k=1,3
        rr=0.0
        do i=1,3
          rr=rr+din(i)*rot(k,i)
        end do
        d(k)=rr
      end do
      call trnsfrd(dout,d,3)
      return
      end
