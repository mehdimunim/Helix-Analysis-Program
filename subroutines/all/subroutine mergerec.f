      subroutine mergerec(line,index,indexo,isegno,iresncol1,iresncol2,
     -  m1,m2,n1,n2,iresn,ireso,isegn,n,maxrec)
c#    MMC routine 127 lstmod: 01/08/88
c*****Merge two sets for the sorting
      character* 132 line(maxrec)
      dimension index(n),indexo(n),isegno(n),iresn(n),ireso(n),isegn(n)
c     i,j: counters for the first and second part, resp.
c     k: counter for the merged array
c     print *,'m1,m2,n1,n2,n=',m1,m2,n1,n2,n
      i=m1
      j=n1
      k=1
20    if (i .le. m2 .and. j .le. n2) then
c       Pick value from first or second part, depending on the comparison
        if (isegno(i) .gt. isegno(j)) then
          idiff=1
        else if (isegno(i) .lt. isegno(j)) then
          idiff=-1
        else
c         Compare residue numbers since segment id's were identical
          call readint(line(index(i)),iresncol1,iresncol2,numi,2,1,
     -      irerr)
          call readint(line(index(j)),iresncol1,iresncol2,numj,2,1,
     -      irerr)
          if (numi .gt. numj) then
            idiff=1
          else
            idiff=0
          end if
        end if
        if (idiff .eq. 1) then
          iresn(k)=index(j)
          ireso(k)=indexo(j)
          isegn(k)=isegno(j)
          j=j+1
        else
          iresn(k)=index(i)
          ireso(k)=indexo(i)
          isegn(k)=isegno(i)
          i=i+1
        end if
        k=k+1
        go to 20
      else
        if (i .le. m2) then
c         Leftover in the first part
          do l=i,m2
            iresn(k)=index(l)
            ireso(k)=indexo(l)
            isegn(k)=isegno(l)
            k=k+1
          end do
        end if
        if (j .le. n2) then
c         Lefotver in the second part
          do l=j,n2
            iresn(k)=index(l)
            ireso(k)=indexo(l)
            isegn(k)=isegno(l)
            k=k+1
          end do
        end if
      end if
c     Move back merged data into index
      k=0
      do l=m1,n2
        k=k+1
        index(l)=iresn(k)
        indexo(l)=ireso(k)
        isegno(l)=isegn(k)
      end do
      return
      end
