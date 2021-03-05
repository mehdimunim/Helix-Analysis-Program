      subroutine restoix(iresix,ir1,ir2,ix1,ix2,ifail,maxres)
      dimension iresix(maxres)
      ifail=0
      ix1=iresix(ir1)
      if (ix1 .eq. 0) then
        ir=ir1
        do while (iresix(ir) .eq. 0 .and. ir .lt. ir2)
          ir=ir+1
        end do
        if (iresix(ir) .gt. 0) then
          ix1=iresix(ir)
        else
          print *,'Residue range [',ir1,',',ir2,'] has no data'
          ifail=1
        end if
      end if
      ix2=iresix(ir2)
      if (ix2 .eq. 0) then
        ir=ir2
        do while (iresix(ir) .eq. 0 .and. ir .gt. ir1)
          ir=ir-1
        end do
        ix2=iresix(ir)
        if (iresix(ir) .eq. 0)print *,'PROGRAM ERROR in restoix'
      end if
c     print *,'RESTOIX ir1,ir2=',ir1,ir2,' ix1,ix2=',ix1,ix2
      return
      end
