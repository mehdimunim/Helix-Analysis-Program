      subroutine bonddescr(ib,ihbpair,line,index,iresno,isegno,
     -  inamcol1,inamcol2,irescol1,irescol2,bond,lbond,bondlab,
     -  lbondlab,nbb,ianc_anc,isc,ia1,ia2,ir1,ir2,maxrec,maxbonds)
      dimension ihbpair(2,maxbonds),index(maxrec),iresno(maxrec),
     -  isegno(maxrec),ianc_anc(maxbonds),isc(maxrec)
      character*132 line(maxrec)
      character*(*) bondlab
      character*2 sc_lab(2,2)
      character*(*) bond
      data sc_lab /'BB','SB','BS','SS'/
      ia1=ihbpair(1,ib)
      ia2=ihbpair(2,ib)
      ir1=iresno(ia1)
      ir2=iresno(ia2)
      write (bond,1000)
     -    line(index(ia1))(inamcol1:inamcol2),ia1,
     -    line(index(ia1))(irescol1:irescol2),ir1,isegno(ia1),
     -    line(index(ia2))(inamcol1:inamcol2),ia2,
     -    line(index(ia2))(irescol1:irescol2),ir2,isegno(ia2)
      lbond=2*(inamcol2-inamcol1+irescol2-irescol1+2)+38
      lbondlab=1
      bondlab=' '
      if (nbb .gt. 0) then
        bondlab(lbondlab+1:lbondlab+3)=
     -    sc_lab(isc(ihbpair(1,ib))+1,isc(ihbpair(2,ib))+1)
          lbondlab=lbondlab+2
        end if
        if (ianc_anc(ib) .eq. 1) then
          bondlab(lbondlab+1:lbondlab+3)=' AA'
          lbondlab=lbondlab+3
        end if
      return
1000  format(1x,a,i6,1x,a,i4,' C/S',i2,' - ',a,i6,1x,a,i4,' C/S',i2)
      end
