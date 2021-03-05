      subroutine findCA(line,index,ifr,ilr,inamcol1,inamcol2,iCA,nnoCA,
     -  maxrec)
      dimension index(maxrec)
      character* 132 line(maxrec)
      character*8 aname
      iCA=0
      do ia=ifr,ilr
        aname(1:inamcol2-inamcol1+1)=line(index(ia))(inamcol1:inamcol2)
        if (aname(1:4) .eq. 'CA  ' .or. aname(1:4) .eq. ' CA ') iCA=ia
      end do
      if (iCA .eq. 0) then
        iCA=ifr
        nnoCA=nnoCA+1
      end if
      return
      end
