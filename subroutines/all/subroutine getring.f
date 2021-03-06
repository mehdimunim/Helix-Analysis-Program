      subroutine getring(line,index,ix5,irescol1,irescol2,
     -  inamcol1,inamcol2,n,numres,iresring,iresno,ifres,ilres,nmem,
     -  ansrun,noask,incgen,iapex,maxring,maxrec)
      character* 132 line(maxrec)
      dimension index(maxrec),ix5(maxring),iresno(n),ifres(maxrec),
     -  ilres(maxrec)
      character*1 ansrun
      character*2 s5name(5),pname(5),s6name(6)
      character*4 atomnam
      character*15 question
      data s5name /'C1','C2','C3','C4','O4'/
      data s6name /'C1','C2','C3','C4','C5','O5'/
      data pname /'CA','CB','CG','CD','N '/
      data question /'Ring member #  '/
      notfound=0
100   if (noask .eq. 0)
     -   call quiz(ansrun,irtyp,' ',' ',0,'ring type',9,0,5,6,0)
      if (ansrun .eq. 'd') then
        nmem=0
        return
      else if (ansrun .eq. 'g') then
        call getint('Number of atoms in the ring',27,nmem,1,MAXRING,
     -    nmem,0)
        do i=1,nmem
          write (question(14:15),1000) i
          call getint(question,15,999999,1,n,ix5(i),00)
        end do
        incgen=0
      else if (noask .eq. 0) then
        call getint('Residue number of the ring',26,0,1,numres,iresring,
     -    0)
        if (iresring .lt. 1 .or. iresring .gt. iresno(n)) then
          write (6,1001) iresring,iresno(n)
          go to 100
         end if
      end if
      if (ansrun .eq. '6') then
c       Find hexose sugar
        nmem=6
c       Apex is O
        incgen=5
        do i=1,nmem
          do ia=ifres(iresring),ilres(iresring)
            atomnam=line(index(ia))(inamcol1:inamcol1+3)
            call leftadjust4(atomnam,atomnam)
            if (atomnam(1:2) .eq. s6name(i) .and.
     -          atomnam(3:3) .ne. ' ') then
              ix5(i)=ia
              go to 110
            end if
          end do
          write (6,1002) s6name(i),iresring,
     -      line(index(ifres(iresring)))(irescol1:irescol2)
          notfound=1
110       continue
        end do
      else if (ansrun .eq. '5') then
c       Find pentose sugar
        nmem=5
c       Apex is O
        incgen=4
        do i=1,nmem
          do ia=ifres(iresring),ilres(iresring)
            atomnam=line(index(ia))(inamcol1:inamcol1+3)
            call leftadjust4(atomnam,atomnam)
            if (atomnam(1:2) .eq. s5name(i) .and.
     -          atomnam(3:3) .ne. ' ') then
              ix5(i)=ia
              go to 120
            end if
          end do
          notfound=1
          write (6,1002) s5name(i),iresring,
     -      line(index(ifres(iresring)))(irescol1:irescol2)
120       continue
        end do
      else if (ansrun .eq. 'p') then
c       Find proline
        nmem=5
c       Apex is N (?)
        incgen=4
        do i=1,nmem
          do ia=ifres(iresring),ilres(iresring)
c           print *,'ia,r=',ia,line(index(ia))(inamcol1:inamcol1+1)
            atomnam=line(index(ia))(inamcol1:inamcol1+3)
            call leftadjust4(atomnam,atomnam)
            if (atomnam(1:2) .eq. pname(i)) then
              ix5(i)=ia
              go to 130
            end if
          end do
          notfound=1
          write (6,1002) pname(i),iresring,
     -      line(index(ifres(iresring)))(irescol1:irescol2)
130       continue
        end do
      end if
      nnamcol=inamcol2-inamcol1+1
      if (notfound .eq. 1) then
        nmem=0
      else
        if (nnamcol .gt. 0) then
        write (6,1003) (ix5(i),line(index(ix5(i)))(irescol1:irescol2),
     -    line(index(ix5(i)))(inamcol1:inamcol2),i=1,nmem)
        iapex=ix5(incgen+1)
        write (6,1007) iapex,' ',
     -    line(index(ix5(incgen+1)))(irescol1:irescol2),
     -    line(index(ix5(incgen+1)))(inamcol1:inamcol2)
        else
          write (6,1004) (ix5(i),i=1,nmem)
          write (6,1007) iapex
        end if
      end if
      return
1000  format(i2)
1001  format(' ERROR: residue number ',i5,' is outside the [0,',i5,
     -  '] range')
1002  format(' ERROR: ring atom ',a,' is not found in residue ',i4,
     -  ' (',a,')')
1003  format(' Calculations for ring',/,5(i6,' ',a,1x,a))
1004  format(' Calculations for ring',/,(5i6))
1007  format(' Apex of ring: atom',i6,a,'(',a,1x,a,')')
      end
