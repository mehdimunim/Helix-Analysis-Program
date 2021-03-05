      subroutine fitpoints(c,n,ndim,c0,dir,axfact,idebughx)
      real*8 c(3,n),c0(3),dir(3),axfact(n)
      real*8 dd(3),sum,ddot
      real*8 rr(3,3),diag(3),offdiag(3)
      data kcol /0/,kmin /0/,kmax /0/
c     For n points in c, find the mean (c0) and the direction of the normal
c     to the best fitting plane (ndim=3) or the direction of the best fitting
c     line (ndim=2)
c     print *,'FITPOINTS n=',n
      if (ndim .ne. 2 .and. ndim .ne. 3) then
        print *,'PROGRAM ERROR: invalid ndim in fitpoints=',ndim
        stop
      end if
      call zeroitd(c0,3)
      do i=1,n
        do k=1,3
          c0(k)=c0(k)+c(k,i)
        end do
      end do
      do k=1,3
        c0(k)=c0(k)/dfloat(n)
      end do
      if (idebughx .eq. 2) return
      do k=1,3
        do l=k,3
          sum=0.d0
          do i=1,n
            sum=sum+(c(k,i)-c0(k))*(c(l,i)-c0(l))
          end do
          rr(k,l)=sum
          rr(l,k)=sum
        end do
      end do
      if (idebughx .gt. 0) write (77,1000) 'Matrix',rr
c     Find the eigenvectors a and eigenvalues mu of rr
      call dtred2(rr,3,3,diag,offdiag)
      call dtqli(diag,offdiag,3,3,rr,ierr)
      if (ierr .gt. 0) then
        write (6,1004)
        stop
      end if
      if (idebughx .gt. 0) write (77,1000) 'Eigenvectors',rr
      if (idebughx .gt. 0) then
        write (77,1000) 'diag',diag
        write (77,1000) 'offdiag',offdiag
      end if
c     The columns of the matrix rr are the eigenvectors
c     The eigenvalues are in diag
      emax=0.0
      emin=1000.0
      nz=0
      do k=1,3
        if (dabs(diag(k)) .lt. 0.005d0) then
          nz=nz+1
        else
          if (diag(k) .gt. emax) then
            emax=diag(k)
            kmax=k
          end if
          if (diag(k) .lt. emin) then
            emin=diag(k)
            kmin=k
          end if
        end if
      end do
      if (ndim .eq. 2 .and. nz .eq. 0) write (6,1001) diag
      if (ndim .eq. 3 .and. nz .gt. 0) print *,'WARNING: 3D dataset ',
     -  'appears to be in a plane'
      if (ndim .eq. 2) kcol=kmax
      if (ndim .eq. 3) kcol=kmin
      do k=1,3
        dir(k)=rr(k,kcol)
      end do
      do i=1,n
        call dvdif(c(1,i),c0,dd)
        axfact(i)=ddot(dir,dd)
      end do
      return
1000  format(1x,a,'=',/,(3f15.5))
1001  format(' ERROR: no zero eigenvalue found when fitting points ',
     -  'in a plane',/,8x,'diag=',3f10.5)
1004  format(' Calculation aborted due to diagonalization failure')
      end
