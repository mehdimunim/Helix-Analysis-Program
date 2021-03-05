      subroutine condensemask(index,ires1,ires2,iout,mxrsd)
      dimension index(mxrsd)
      character*80 line
c     print *,'CONDENSEMASK ires1,ires2=',ires1,ires2,' iout=',iout
c      write (6,8888) (index(i),i=1,n)
c8888  format(' index:',/,(10i5))
      ic=1
      i=ires1
      do while (i .le. ires2)
c       Find first/next nonzero
        do while (index(i) .eq. 0 .and. i .le. ires2)
          i=i+1
        end do
        if (i .le. ires2) then
c         Find out if singleton or range
          i1=i
          do while (index(i) .gt. 0 .and. i .le. ires2)
            i=i+1
          end do
          i2=i-1
c         print *,'i1=',i1,' i2=',i2,' ic=',ic
          call printrange(line,i1,i2,ic,0,iout)
        end if
      end do
      if (ic .gt. 1) write (iout,1000) line(1:ic-2)
      return
1000  format(1x,a)
      end
