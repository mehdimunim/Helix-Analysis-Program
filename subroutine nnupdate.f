      subroutine nnupdate(index,indexn,indexo,n,nmax,line,isort,maxrec)
      dimension index(nmax),indexn(nmax),indexo(nmax)
      character* 132 line(maxrec)
c     Change atomnumbers in bond list
      do ia=1,n
        do ib=1,6
          ic1=6+(ib-1)*8
          ic2=ic1+4
          read (line(index(ia))(ic1:ic2),1000) iold
c         if (iold .ne. 0) write (6,888) ia,index(ia),ib,ic1,ic2,iold
c888      format(' ia,index(ia),ib=',3i3,' ic1,2=',2i3,' iold=',i3)
          if (iold .gt. 0) then
c           print *,'isort,indexo(iold),indexn(iold)=',
c    -         isort,indexo(iold),indexn(iold)
            if (isort .eq. 0) then
              if (iold .gt. 0)
     -          write (line(index(ia))(ic1:ic2),1000) indexn(iold)
            else
              if (indexo(iold) .gt. 0)
     -          write (line(index(ia))(ic1:ic2),1000) indexo(iold)
            end if
          end if
        end do
      end do
      return
1000  format(i5)
      end
