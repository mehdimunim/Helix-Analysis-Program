      subroutine checkforhelix(nhxres,dssplab,indexn,iw,hxoklab,ihxok,
     -  lab,llab,maxrsd,maxrec)
      dimension indexn(maxrec)
      character*1 dssplab(maxrsd)
      character*6 hxoklab(3)
      character*(*) lab
      ihxok=1
      do ir=1,nhxres
        indexn(ir)=1
        if (dssplab(ir) .ne. 'G' .and. dssplab(ir) .ne. 'H'
     -    .and. dssplab(ir) .ne. 'I' .and.
     -     dssplab(ir) .ne. 'L') indexn(ir)=0
        if (dssplab(ir) .eq. ' ') dssplab(ir)='?'
      end do
      if (indexn(1) .eq. 0 .or. indexn(nhxres) .eq. 0) then
        ihxok=3
      else
        do ir=2,nhxres-1
          if (indexn(ir) .eq. 0) ihxok=2
        end do
      end if
      write (iw,2088) lab(1:llab),hxoklab(ihxok),
     -  (dssplab(ir),ir=1,nhxres)
      return
2088  format(a,' is ',a,': ',(60a1))
      end
