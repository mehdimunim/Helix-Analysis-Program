      subroutine setcol(inpcrdtyp,ncol,idcol,ialtcol,iinscol,
     -  inamcol1,inamcol2,irescol1,irescol2,iccol1,iccol2,
     -  iresncol1,iresncol2,iseqncol1,iseqncol2,isegcol1,isegcol2,
     -  iresidcol1,iresidcol2,iqcol1,iqcol2,ipotcol1,ipotcol2,
     -  iocccol1,iocccol2,ichemcol1,ichemcol2,nrescol,nresncol,nsegcol,
     -  nnamcol,iofull)
      common /columnlim/ incol(19),iidcol(19),iialtcol(19),iiinscol(19),
     -  iinamcol(2,19),iirescol(2,19),iiccol(2,19),iiresncol(2,19),
     -  iiseqncol(2,19),iisegcol(2,19),iiresidcol(2,19),iiqcol(2,19),
     -  iipotcol(2,19),iiocccol(2,19),iichemcol(2,19)
c     print *,'SETCOL inpcrdtyp=',inpcrdtyp
      ncol=incol(inpcrdtyp)
      idcol=iidcol(inpcrdtyp)
      ialtcol=iialtcol(inpcrdtyp)
      iinscol=iiinscol(inpcrdtyp)
      inamcol1=iinamcol(1,inpcrdtyp)
      inamcol2=iinamcol(2,inpcrdtyp)
      irescol1=iirescol(1,inpcrdtyp)
      irescol2=iirescol(2,inpcrdtyp)
      iccol1=iiccol(1,inpcrdtyp)
      iccol2=iiccol(2,inpcrdtyp)
      iresncol1=iiresncol(1,inpcrdtyp)
      iresncol2=iiresncol(2,inpcrdtyp)
      iseqncol1=iiseqncol(1,inpcrdtyp)
      iseqncol2=iiseqncol(2,inpcrdtyp)
      isegcol1=iisegcol(1,inpcrdtyp)
      isegcol2=iisegcol(2,inpcrdtyp)
      iresidcol1=iiresidcol(1,inpcrdtyp)
      iresidcol2=iiresidcol(2,inpcrdtyp)
      iqcol1=iiqcol(1,inpcrdtyp)
      iqcol2=iiqcol(2,inpcrdtyp)
      ipotcol1=iipotcol(1,inpcrdtyp)
      ipotcol2=iipotcol(2,inpcrdtyp)
      iocccol1=iiocccol(1,inpcrdtyp)
      iocccol2=iiocccol(2,inpcrdtyp)
      ichemcol1=iichemcol(1,inpcrdtyp)
      ichemcol2=iichemcol(2,inpcrdtyp)
      nrescol=irescol2-irescol1+1
      nresncol=iresncol2-iresncol1+1
      nsegcol=isegcol2-isegcol1+1
      nnamcol=inamcol2-inamcol1+1
      if (inpcrdtyp .gt. iofull) then
        nrescol=4
        iqcol1=0
      end if
      return
      end
