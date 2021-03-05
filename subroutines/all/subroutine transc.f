      subroutine transc(c0,del,rlambda)
      dimension c0(3),del(3)
      do k=1,3
        c0(k)=c0(k)+rlambda*del(k)
      end do
      return
      end
