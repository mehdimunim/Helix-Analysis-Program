      subroutine rounddiv(max,maxdiv,ndiv,idiv)
      dimension limits (27)
      data limits /1,2,5,24*1/
      do i=1,8
        do k=1,3
          limits(3*i+k)=10*limits(3*(i-1)+k)
        end do
      end do
      do i=1,27
        if (max .le. limits(i)*maxdiv) then
          idiv=limits(i)
          ndiv=(max-1)/idiv +1
          return
        end if
      end do
      idiv=limits(27)
      ndiv=(max-1)/idiv +1
      return
      end
