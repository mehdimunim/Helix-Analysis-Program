      subroutine unitvec(r1,r2,e)
      dimension r1(3),r2(3),e(3)
      call arrdiff(r2,r1,e,3)
      call norm(e,1.0)
      return
      end
