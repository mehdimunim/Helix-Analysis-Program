      subroutine scramble(ix,n)
      dimension ix(n)
      dimension rand(1)
      call indexit(ix,1,n,0)
      left=n
      do i=1,n
        call randpx(1,rand)
        is1=i+rand(1)*left
        ix0=ix(is1)
        ix(is1)=ix(i)
        ix(i)=ix0
        left=left-1
      end do
      return
      end
