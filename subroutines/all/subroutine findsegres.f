      subroutine findsegres(isegno,iresno,ixres,ifst,ilst,iseg,ires,
     -  iaf,irf,ifail)
c*****Find the residue sequence number of segid=iseg, resid=ires
      dimension isegno(ilst),iresno(ilst),ixres(ilst)
      ifail=0
c     print *,'FINDSEGRES ifst,ilst=',ifst,ilst
      do ia=ifst,ilst
        if (isegno(ia) .eq. iseg) then
          if (iresno(ia) .eq. ires) then
            iaf=ia
            irf=ixres(ia)
c           print *,'FINDSEGRES iseg,ires,ia,irf=',iseg,ires,ia,irf
            return
          end if
        end if
      end do
      write (6,1000) iseg,ires,ifst,ilst
      ifail=1
      return
1000  format(' ERROR: segment number',i3,' residue ',i4,' was not ',
     -  'found',/,' Atom range searched:',i5,' - ',i5)
      end
