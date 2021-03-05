      subroutine mergelst(indexx,value,m1,m2,n1,n2,ires,res,maxt)
c#    MMC routine 263 lstmod: 06/23/97
c*****Merge two sets for the sorting
      dimension indexx(maxt),value(maxt),ires(maxt),res(maxt)
c     i,j: counters for the first and second part, resp.
c     k: counter for the merged array
      i=m1
      j=n1
      k=1
20    if (i .le. m2 .and. j .le. n2) then
c       Pick value from first or second part, depending on value
        if (value(i) .le. value (j)) then
          res(k)=value(i)
          ires(k)=indexx(i)
          i=i+1
          go to 40
        end if
        res(k)=value(j)
        ires(k)=indexx(j)
        j=j+1
40      k=k+1
        go to 20
      end if
      if (i .le. m2) then
c       Leftover in the first part
        do l=i,m2
          res(k)=value(l)
          ires(k)=indexx(l)
          k=k+1
        end do
      end if
      if (j .le. n2) then
c       Leftover in the second part
        do l=j,n2
          res(k)=value(l)
          ires(k)=indexx(l)
          k=k+1
        end do
      end if
c     Move back merged data into index, value
      k=0
      do l=m1,n2
        k=k+1
        indexx(l)=ires(k)
        value(l)=res(k)
      end do
      return
      end
