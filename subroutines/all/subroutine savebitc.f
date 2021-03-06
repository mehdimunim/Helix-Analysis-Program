      subroutine savebitc(mapbit,ibitx,n,nbits,maxbit)
c*****Save the bits from ibitx into mapbit
      dimension mapbit(maxbit),ibitx(n)
      iw=1
      ib=0
      ibb=0
      mapbi=0
      do i=1,n
        ib=ib+1
        ibb=ibb+1
        if (ibitx(ibb) .eq. 1) mapbi=ibset(mapbi,ib-1)
        if (ib .eq. nbits) then
c       write (77,*)' write iw=',iw,' mapbi=',mapbi,
c    -    ' ibdone,itodo=',ibdone,itodo
          mapbit(iw)=mapbi
          iw=iw+1
          mapbi=0
          ib=0
        end if
      end do
      if (ib .gt. 0) mapbit(iw)=mapbi
c     if (ib .gt. 0) write (77,*)' writ iw=',iw,' mapbi=',mapbi
      return
      end
