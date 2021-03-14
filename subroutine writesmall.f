      subroutine writesmall(iout,a,n,iy,iyorig)
      dimension a(n)
      character*80 line
c     print *,'WRITESMALL n,iy,iyorig=',n,iy,iyorig
      nw=0
      do while (nw .lt. n)
        call blankout(line,1,80)
        if (nw .eq. 0) write (line(1:24),1001) iy,iyorig
        write (line(25:74),1002) (a(i),i=nw+1,min0(n,nw+10))
        do i=1,min0(10,n-nw)
          if (a(nw+i) .lt. 0.00001) line(24+(i-1)*5+1:24+i*5)=' 0.0 '
        end do
        write (iout,1000) line(1:74)
        nw=nw+10
      end do
      return
1000  format(a)
1001  format(' iy=',i4,' iy(orig)=',i5,':')
1002  format(10f5.2)
      end
