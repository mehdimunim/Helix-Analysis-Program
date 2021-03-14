      subroutine findat(iarep,ira1,ira2,line,index,irescol1,inamcol1,
     -  maxrec)
c*****Find the representative atom of residue ir
      character* 132 line(maxrec)
      dimension index(maxrec)
      character*3 resnam3,repats,repat0
      character*4 atname,an
      character*8 resnam,rn
      common /represent/ repats(2,50),maxrepat
c     Find first atom of residue ir
c     write (77,*)'FINDAT START ir,ira1,ira2=',ir,ira1,ira2
c     Find residue name and corresponding representative atom
      resnam='        '
      resnam=line(index(ira1))(irescol1:irescol1+4)
      call leftadjustn(resnam,rn,8)
      resnam3=rn(1:3)
      infound=0
c     print *,'resnam3=',resnam3,' resnam=',resnam
      do i=1,maxrepat
        if (resnam3 .eq. repats(1,i)) infound=i
      end do
c     Find representative atoms of residue ir
      iarep=0
      repat0='CA '
      do ia=ira1,ira2
        if (infound .gt. 0) repat0=repats(2,infound)
        atname=line(index(ia))(inamcol1:inamcol1+3)
        call leftadjust4(atname,an)
        if (an(1:3) .eq. repat0) iarep=ia
      end do
      if (iarep .eq. 0 .and. infound .gt. 0)
     -  write (6,1000) repats(2,infound),repats(1,infound),ira1,ira2
c     if (iarep .eq. 0 .and. infound .gt. 0) then
c     do ia=ira1,ira2
c       print *,'ia=',ia
c       print *,line(index(ia))(1:78)
c     end do
c     stop
c     end if
      if (iarep .eq. 0) iarep=ira1
      return
1000  format(' Representative atom ',a,' for residue ',a,' not found ',
     - 'in atom range',i6,' - ',i6)
      end
