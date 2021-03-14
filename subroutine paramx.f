      subroutine paramx(cent,ax,a,x)
      dimension cent(3),ax(3),x(3)
      do k=1,3
       x(k)=cent(k)+a*ax(k)
      end do
      return
      end
