      subroutine get_heavyat(ia,nneig,ineig,ixres,nframe,ihb0,maxneig,
     -  maxrec)
      dimension nneig(maxrec),ineig(maxneig,maxrec),ixres(maxrec)
c     print *,'GET_HEAVYAT ia=',ia,' MAXNEIG,MAXREC=',maxneig,maxrec
      ihb0=0
      nx=0
      ing_prev=0
      do iaa=1,nneig(ia)
        ing=ineig(iaa,ia)
        if (nneig(ia) .eq. 1 .or. ixres(ia) .eq. ixres(ing)) then
          ihb0=ineig(iaa,ia)
          ing_prev=ing
          nx=nx+1
          if (nx .gt. 1) then
            if (nx .eq. 2)
     -        write (6,2002) ia,ixres(ia),ing,ixres(ing_prev)
            write (6,2002) ia,ixres(ia),ing,ixres(ing)
          end if
        end if
      end do
      if (nx .gt. 1) then
        if (nframe .eq. 0) write (6,2001) ia,nx
        if (nframe .gt. 0) write (6,2001) ia,nx,' ',nframe
      end if
      return
2001  format(' WARNING: Hbond donor H ',i6,' has ',i3,' bond(s)',a,
     -  ' Nframe=',i6)
2002  format(' H atom',i6,' ixres',i5,' is bonded to atom',i6,' ixres',
     -  i5) 
      end
