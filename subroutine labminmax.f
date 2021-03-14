      subroutine labminmax(ir,ncol,ixcolmin,ixcolmax,minmaxlab1,mxrsd)
      dimension ixcolmin(mxrsd),ixcolmax(mxrsd)
      character*1 minmaxlab1(mxrsd)
c     print *,'LABMINMAX ir,ncol=',ir,ncol
      do ic=1,ncol
        minmaxlab1(ic)=' '
        if (ixcolmax(ic) .eq. ir) then
          minmaxlab1(ic)='M'
        end if
        if (ixcolmin(ic) .eq. ir) then
          minmaxlab1(ic)='m'
        end if
        if (ixcolmin(ic) .eq. ixcolmax(ic) .and. ixcolmax(ic) .eq. ir)
     -    minmaxlab1(ic)='c'
      end do
      return
      end
