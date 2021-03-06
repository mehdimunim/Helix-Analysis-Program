      subroutine writeitem(line,icol1,icol2,name,namlen)
c     Place one item into a record
      character* 132 line
      character*(*) name
      if (icol1 .gt. icol2) return
      if (namlen .eq. 0) then
        namelen=icol2-icol1+1
      else
        namelen=min0(icol2-icol1+1,namlen)
      end if
      call blankout(line,icol1,icol2)
      line(icol1:icol1+namelen-1)=name(1:namelen)
      return
      end
