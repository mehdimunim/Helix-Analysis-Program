      subroutine extension(c,ih,nnh,nstart,n,cmin,cmax,c0,iprint,
     -  ivcheck,vol)
      dimension c(3,n),ih(n),cmin(3),cmax(3),c0(3)
      character*1 xyz
      common /axislab/ xyz(3)
c     Find an approximate center first
c     write (06,*) 'EXTENSION nstart,n=',nstart,n,' nnh=',nnh
      do k=1,3
        cmin(k)=1.e+25
        cmax(k)=-1.e+25
      end do
      do k=1,3
        if (nnh .gt. 0) then
          do ii=nstart,nnh
            i=ih(ii)
            if (c(k,i) .lt. cmin(k)) cmin(k)=c(k,i)
            if (c(k,i) .gt. cmax(k)) cmax(k)=c(k,i)
c           if (iprint .gt. 1) write (77,7711) ii,i,k,c(k,i)
          end do
        else
          do i=nstart,n
            if (c(k,i) .lt. cmin(k)) cmin(k)=c(k,i)
            if (c(k,i) .gt. cmax(k)) cmax(k)=c(k,i)
c           if (iprint .gt. 1) write (77,7711) i,i,k,c(k,i)
          end do
        end if
        c0(k)=(cmax(k)+cmin(k))/2.0
c7711    format(2i5,' c',i1,'=',f10.5)
      end do
      vol=1.0
      do k=1,3
        vol=vol*(cmax(k)-cmin(k))
      end do
      if (iprint .gt. 0)
     -  write (6,1000) (xyz(k),cmin(k),c0(k),cmax(k),k=1,3),vol
      if (ivcheck .eq. 1) then
        do k=1,3
          imax=cmax(k)
          if (imax .eq. 9999) write (6,1001) xyz(k)
        end do
        if (vol/n .gt. 1000.0) then
          write (6,1002) vol,n
          if (iprint .eq. 0)
     -      write (6,1000) (xyz(k),cmin(k),c0(k),cmax(k),k=1,3)
          call askstop(0)
        end if
      end if
      return
1000  format((' Smallest, middle and largest ',a1,
     -  ' coordinate values=',3f9.4,/),
     -  ' Volume of enclosing rectangle=',f12.2,' A^3')
1001  format(' WARNING: input structure contains ',a1,' coordinate ',
     -  'value(s) of 9999.0',/,10x,'indicating undefined values')
1002  format(' WARNING: initial enclosing rectangle volume (',e12.4,
     -  ' for',i7,' atoms)',/,8x,'suggests corrupted input structure')
      end
