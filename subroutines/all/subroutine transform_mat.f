      subroutine transform_mat(r,n,nmax,rijrmin,iuout)
      dimension r(nmax,nmax)
      character*1 ans
      common /transform_dist_dat/ itranstyp,nexp,rijmin,fracexp
      fracexp=1.0
      rijmin=rijrmin
      call quiz(ans,itranstyp,'u',' ',0,'transformation type',19,
     -   0,5,6,00)
      if (itranstyp .gt. 1) then
        if (itranstyp .ne. 4) then
          call getint('Exponent',8,1,1,6,nexp,0)
          fracexp=1.0/float(nexp)
        end if
        do i=1,n
          do j=1,i-1
            call transform_dist(r(i,j),r(i,j))
            r(j,i)=r(i,j)
          end do
        end do
        write (iuout,1001)
        if (itranstyp .eq. 2) write (iuout,1002) nexp
        if (itranstyp .eq. 3) write (iuout,1003) nexp
        if (itranstyp .eq. 4) write (iuout,1004)
        if (itranstyp .eq. 5) write (iuout,1005) nexp
        if (itranstyp .eq. 6) write (iuout,1006) nexp
        if (itranstyp .eq. 7) write (iuout,1007) nexp
        if (itranstyp .eq. 8) write (iuout,1008) nexp
        if (itranstyp .eq. 9) write (iuout,1009) nexp
      else
        write (iuout,*) 'Input matrix used as is'
      end if
      return
1001  format(' Input data trasnformed as')
1002  format(' r(i,j)=r(i,j)**',i1)
1003  format(' r(i,j)=r(i,j)**(1/',i1,')')
1004  format(' r(i,j)=1.0-r(i,j)')
1005  format(' r(i,j)=(1.0-r(i,j))**',i1)
1006  format(' r(i,j)=(1.0-r(i,j))**(1/',i1,')')
1007  format(' r(i,j)=abs(r(i,j)**',i1)
1008  format(' r(i,j)=(1.0-abs(r(i,j))**',i1)
1009  format(' r(i,j)=(r(i,j)-min(r(i,j))**',i1)
      end
