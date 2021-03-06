      subroutine findresnum(iresno,ixres,ir,ifirst,ilast,iafound,
     -  irfound)
      dimension iresno(ilast),ixres(ilast)
c     Find the residue sequence number of residue number ir
c     print *,'if,il=',ifirst,ilast
      ia=ifirst
      do while (iresno(ia) .ne. ir .and. ia .le. ilast)
        ia=ia+1
      end do
      if (ia .gt. ilast) then
        write (6,1000) ir,ifirst,ilast
        ir=0
      end if
      iafound=ia
      irfound=ixres(ia)
      return
1000  format(' ERROR: residue number',i6,' is not found in atom range',
     - i7,' - ',i7)
      end
