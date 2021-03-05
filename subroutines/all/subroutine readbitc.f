      subroutine readbitc(mapbit,ibitx,n,nbits,maxbit)
c*****Extract the bits from mapbit into ibitx
      dimension mapbit(maxbit),ibitx(n)
      nloops=(n-1)/nbits+1
      do il=1,nloops
c       Set loop limits so that inside loops can run in parallel
        ibdone=(il-1)*nbits
        itodo=ibdone+min0(nbits,n-ibdone)
        mapbi=mapbit(il)
c       write (77,*)' read il=',il,' mapbi=',mapbi,
c    -    ' ibdone,itodo=',ibdone,itodo
        do ib=ibdone+1,itodo
          ix=mapbi/2
          ibitx(ib)=mapbi-2*ix
          mapbi=ix
        end do
      end do
      return
      end
