      subroutine sortlist(iout,list,listlen,itemp,iorig,lab,itestonly,
     -  maxval)
      dimension list(listlen),itemp(maxval),iorig(listlen)
      character*3 lab
c     Sort a list of distinct positive integers of maximum value maxval
c     Use it for long lists; listlen about the same as maxval
      call zeroiti(itemp,0,maxval)
      nerr=0
      do i=1,listlen
        if (itemp(list(i)) .gt. 0) then
          write (6,1000) lab,list(i),listlen,maxval
          if (iout .gt. 0) write (iout,1000) lab,list(i),listlen,maxval
          nerr=nerr+1
        end if
        itemp(list(i))=i
      end do
      if (nerr .gt. 0) then
        write (6,1001) (list(i),i=1,listlen)
        if (iout .gt. 0) write (iout,1001) (list(i),i=1,listlen)
      end if
      if (itestonly .eq. 0) then
        ll=0
        do i=1,maxval
          if (itemp(i) .gt. 0) then
            ll=ll+1
            list(ll)=i
            iorig(ll)=itemp(i)
          end if
        end do
      end if
      return
1000  format(' PROGRAM ERROR in ',a,' sortlist: duplicate entry ',i6,
     -  /,' List length=',i6,' Maximum value=',i6)
1001  format(' List entries:',/,(10i6))
      end
