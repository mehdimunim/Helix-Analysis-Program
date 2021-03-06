      subroutine checkbend(c,axisdir,camod,axfact,perpvec,x0,
     -  n,nup,ndown,nrun,nnear,axtol,rcirc,circ,rn,idecide,idebughx)
      real*8 c(3,n),axisdir(3),camod(3,n),perpvec(3,n),axfact(n),x0(3),
     -  circ(3),rn(3)
      real*8 rm(3),x(3),x1(3),x2(3),rr,rmin,ddistsq,ddot
      real*8 cx(3,50)
      dimension nupdown(2)
c     print *,'CHECKBEND n=',n
      call fitpoints(c,n,3,x0,rn,axfact,idebughx)
      if (idebughx .gt. 0) then
        do i=1,n
          call writepdbd(77,c(1,i),i,i,'C   ','CA ','A',1.0,0.0)
        end do
        call writepdbd(77,x0,n+1,n+1,'S   ','X0 ','B',1.0,0.0)
        if (idebughx .gt. 1) then
          idecide=6
          nup=0
          ndown=0
          nrun=0
        else
          call dvsum(x0,rn,x1)
          call writepdbd(77,x1,n+2,n+2,'S   ','RN ','B',1.0,0.0)
        end if
      end if
      if (idebughx .gt. 0)
     -  write (78,1000) (i,(c(k,i),k=1,3),(camod(k,i),k=1,3),axfact(i),
     -    (perpvec(k,i),k=1,3),i=1,n)
c     Find the shortest axis-Calpha distance
      rmin=1000.d0
      do i=1,n
        rr=ddistsq(c(1,i),camod(1,i))
        if (rr .lt. rmin) then
          rmin=rr
          imin=i
        end if
      end do
      rr=dsqrt(rmin)
c     "Pull" uniformly all the Calphas toward the axis as close as possible
      do i=1,n
        do k=1,3
          camod(k,i)=c(k,i)+rr*perpvec(k,i)
        end do
      end do
      if (idebughx .gt. 0) then
        do i=1,n
          call writepdbd(77,camod(1,i),i,i,'C   ','CA ','B',1.0,0.0)
        end do
      end if
c     Calculate the radius of the fitting circle to the pulled points
      call circfit(c,n,axisdir,circ)
      rcirc=dsqrt(ddistsq(circ,x0))
      rr=0
      do i=1,n
        rr=rr+ddistsq(circ,c(1,i))
      end do
      rcirc=dsqrt(rr/dfloat(n))
      if (idebughx .gt. 0) write (77,1002) circ,x0,rcirc
c     Find the projection of the pulled points on the fitting plane
      do i=1,n
        call dvdif(camod(1,i),x0,x1)
        rr=ddot(x1,rn)
        do k=1,3
          camod(k,i)=camod(k,i)-rr*rn(k)
        end do
      end do
      call fitpoints(camod,n,2,x0,rm,axfact,0)
      if (idebughx .gt. 0) then
        do i=1,n
          call writepdbd(77,camod(1,i),i,i,'N   ','CAP','C',1.0,0.0)
        end do
        call writepdbd(77,x0,n+1,n+1,'P   ','X0 ','B',1.0,0.0)
        call dvsum(x0,rm,x1)
        call writepdbd(77,x1,n+2,n+2,'P   ','RM ','B',1.0,0.0)
      end if
      if (idebughx .gt. 0) then
c       Check if projections are indeed in a plane
        call dvdif(camod(1,1),x0,x1)
        call dvdif(camod(1,2),x0,x2)
        call dvprd(x1,x2,x)
        do i=3,n
          call dvdif(camod(1,i),x0,x1)
          rr=ddot(x,x1)
          write (78,1001) i,(camod(k,i),k=1,3),rr
        end do
        call trnsfrd(cx,camod,3*n)
      end if
c     Get the distance vectors between the point and its projection on the line
      do i=1,n
        do k=1,3
          camod(k,i)=(camod(k,i)-x0(k))-axfact(i)*rm(k)
        end do
      end do
      if (idebughx .gt. 0) then
        do i=1,n
          call dvdif(cx(1,i),camod(1,i),x)
          call writepdbd(77,x,i,i,'O   ','CA0','D',1.0,0.0)
        end do
      end if
      nnear=0
      if (axtol .gt. 0.0) then
c       Eliminate (and count) residuea that are within tolerance of the axis
        do i=1,n
          if (sngl(ddot(camod(1,i),camod(1,i))) .le. axtol**2) then
            nnear=nnear+1
          else if (nnear .gt. 0) then
            call trnsfrd(camod(1,i-nnear),camod(1,i),3)
          end if
        end do
      end if
      call zeroiti(nupdown,0,2)
      iupdown=0
      nswitch=0
      do i=2,n-nnear
        if (ddot(camod(1,i-1),camod(1,i)) .lt. 0.d0) then
          nswitch=nswitch+1
          iupdown=1-iupdown
        end if
        nupdown(iupdown+1)=nupdown(iupdown+1)+1
      end do
      nrun=nswitch+1
      nup=nupdown(1)
      ndown=nupdown(2)
      call runtest(nup,ndown,nrun,idecide)
      return
1000  format(i3,' c=',3f8.3,' camod=',3f8.3,' axfac=',f8.3,/,
     -  ' perp=',3f8.4)
1001  format(i3,' camod=',3f8.3,' plane check (->0)=',f10.6)
1002  format(' CHECKBEND circ=',3f10.5,' x0=',3f10.5,' rcirc=',f10.5)
      end
