      subroutine mrglimtst(iout,n,maxt,ireturn)
      ireturn=0
      if (n .eq. 1) then
        ireturn=1
      else if (n .gt. maxt) then
        write (iout,1000) n,maxt
        return=1
      else if (n .lt. 0) then
        write (iout,1001) n
        ireturn=1
      end if
      return
1000  format(' ERROR: mrgsrt input array length(',i8,') exceeds the ',
     -  'limit:',i6)
1001  format(' PROGRAM ERROR: negative number of elements to sort:',i6)
c     Individual intervals consist of single elements
      end
