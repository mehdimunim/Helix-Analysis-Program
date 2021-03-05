      subroutine nninit(nneig,nhbneig,ineig,nhneig,nnneig,ncneig,
     -  nsneig,npneig,nfirst,n,ihbinit,maxng)
      dimension nneig(n),ineig(maxng,n),nhbneig(n),nhneig(n),
     -  nnneig(n),ncneig(n),nsneig(n),npneig(n)
      if (ihbinit .eq. 1) then
        do i=nfirst,n
          call zeroiti(ineig(1,i),0,maxng)
          nhbneig(i)=0
        end do
      end if
      do i=nfirst,n
        nhneig(i)=0
        ncneig(i)=0
        nnneig(i)=0
        nsneig(i)=0
        npneig(i)=0
        nneig(i)=0
      end do
      return
      end
