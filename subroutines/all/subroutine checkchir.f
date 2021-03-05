      subroutine checkchir(c,n,i0,i1,i2,i3,i4,isg)
      dimension c(3,n)
      dimension d01(3),d02(3),d34(3),rnorm12(3)
c     Check the chirality of atom c(1,i0). If i1,i2,i3,i4 are the indices of
c     atoms with increasing atomic number then isg = +1 -> R ; isg = -1 -> S.
      call arrdiff(c(1,i1),c(1,i0),d01,3)
      call arrdiff(c(1,i2),c(1,i0),d02,3)
      call arrdiff(c(1,i4),c(1,i3),d34,3)
      call vprd(d01,d02,rnorm12)
      isg=1
      if (scprod(d34,rnorm12) .lt. 0.0) isg=-1
      return
      end
