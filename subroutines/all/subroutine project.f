      subroutine project(x,c,ax,p,a)
      dimension x(3),c(3),ax(3),p(3)
c     Find the nearest point p from x on c+a*ax
      a=scprod(ax,x)-scprod(ax,c)
      do k=1,3
        p(k)=c(k)+a*ax(k)
      end do
      return
      end
