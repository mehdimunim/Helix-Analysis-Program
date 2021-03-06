      subroutine uniquename(anames,n)
      character*4 anames(n)
      dimension ntyp(200)
      character*4 typnames(200),number
      ntyps=0
      do ia=1,n
        call findname(anames(ia),typnames,1,ntyps,ix,4)
        if (ix .eq. 0) then
          ntyps=ntyps+1
          if (ntyps .gt. 200) then
            print *,'Number of names exceeds 200 - redimension ',
     -        'the subroutine uniquenames'
            stop
          end if
          typnames(ntyps)=anames(ia)
        end if
      end do
      call zeroiti(ntyp,0,200)
      do ia=1,n
        call findname(anames(ia),typnames,1,ntyps,ix,4)
        if (ix .eq. 0) then
          print *,'PROGRAM ERROR in uniquename: ix=0'
        else
          ntyp(ix)=ntyp(ix)+1
          call lastchar(anames(ia),lc,4)
          ic=1
          call writeint(number,ic,ntyp(ix),len)
          if (lc+len .gt. 4) then
            write (6,1000) ntyp(ix),anames(ia)
            if (anames(ia)(1:1) .eq. ' ') write (6,1001)
            stop
          end if
          anames(ia)(lc+1:lc+len)=number(1:len)
        end if
      end do
      return
1000  format(' There is no room to write ',i4,' after ',a)
1001  format(' You can make room by first leftadjust the atom names')
      end
