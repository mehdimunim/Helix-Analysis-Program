      subroutine listatom(line,index,iatno,ia,inpcrdtyp,iofull,lab,llab,
     -  n,maxrec)
      dimension index(n),iatno(n)
      character* 132 line(maxrec),ll
      character*(*) lab
      character*4 pflab
      character*12 qlab
      call setcol(inpcrdtyp,ncol,idcol,ialtcol,iinscol,
     -  inamcol1,inamcol2,irescol1,irescol2,iccol1,iccol2,
     -  iresncol1,iresncol2,iseqncol1,iseqncol2,isegcol1,isegcol2,
     -  iresidcol1,iresidcol2,iqcol1,iqcol2,ipotcol1,ipotcol2,
     -  iocccol1,iocccol2,ichemcol1,ichemcol2,nrescol,nresncol,
     -  nsegcol,nnamcol,iofull)
      ll=line(index(ia))
      pflab='    '
      qlab='            '
      nqcol=iqcol2-iqcol1+1
      if (nqcol .gt. 0) qlab(1:nqcol)=ll(iqcol1:iqcol2)
      nqcol=max0(1,nqcol)
      if (ipotcol1 .le. ipotcol2) pflab=ll(ipotcol1:ipotcol2)
      write (6,1000) lab(1:llab),ia,ll(inamcol1:inamcol2),
     -  ll(irescol1:irescol2),ll(iresncol1:iresncol2),
     -  iatno(ia),pflab,qlab(1:nqcol)
      return
1000  format(1x,a,i6,' (',a,1x,a,') resnum=',a,' iatnum=',i2,
     -  ' pf=',a,' q=',a)
      end
