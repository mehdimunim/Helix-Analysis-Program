      subroutine fixrecform(line,index,n,iconv,inamcol1,inamcol2,
     -  irescol1,irescol2,maxrec)
      character* 132 line(maxrec)
      character*4 atnamo,atnamn
      character*8 resnam
      dimension index(n)
      if (iconv .eq. 2) itofrom=1
c     'Regularize' PDB atomnames
      if (iconv .eq. 3) itofrom=0
c     'De-regularize' PDB atomnames
      n2=0
      do ia=1,n
        resnam='        '
        resnam(1:irescol2-irescol1+1)=
     -    line(index(ia))(irescol1:irescol2)
        atnamo=line(index(ia))(inamcol1:inamcol2)
        call regularpdb(atnamo,atnamn,itofrom)
        line(index(ia))(inamcol1:inamcol2)=atnamn
      end do
      if (n2 .gt. 0) print *,'WARNING: ',n2,
     -    ' atomnames had more than one leading digits'
      return
      end
