      subroutine batchmean(npts,nskip,u,label,llabel,iout,nopr,avg,sd,
     -  ci)
c#    MMC routine 255 lstmod: 05/02/13
c*****Computes error bound with the method of batch means
      character*(*) label
      dimension u(npts)
      parameter (MAXSORT=128)
      common /sortsat/ ixdat(MAXSORT),dat(MAXSORT),datsort(MAXSORT),
     -  ifst(MAXSORT),ilst(MAXSORT),itemp(MAXSORT),temp(MAXSORT)
      real*8 datsum,datsum2
      dimension nmncrt(20,20),nmxcrt(20,20)
      character*12 uncorr,corr,low,decide
      data uncorr/'Uncorrelated'/,corr/'Correlated  '/,
     -  low/' ???        '/
c     Minimum critical values for correlation test:
      data nmncrt/20*0,
     -  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2,
     -  0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3,
     -  0, 0, 0, 0, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4,
     -  0, 0, 0, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5,
     -  0, 0, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6,
     -  0, 0, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6,
     -  0, 0, 2, 3, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7,
     -  0, 0, 2, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8,
     -  0, 0, 2, 3, 3, 4, 5, 5, 5, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9,
     -  0, 0, 2, 3, 4, 4, 5, 5, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 9,
     -  0, 2, 2, 3, 4, 4, 5, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9,10,10,
     -  0, 2, 2, 3, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 9,10,10,10,10,
     -  0, 2, 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 9, 9,10,10,10,11,11,
     -  0, 2, 3, 3, 4, 5, 6, 6, 7, 7, 8, 8, 9, 9,10,10,11,11,11,12,
     -  0, 2, 3, 4, 4, 5, 6, 6, 7, 8, 8, 9, 9,10,10,11,11,11,12,12,
     -  0, 2, 3, 4, 4, 5, 6, 7, 7, 8, 9, 9,10,10,11,11,11,12,12,13,
     -  0, 2, 3, 4, 5, 5, 6, 7, 8, 8, 9, 9,10,10,11,11,12,12,13,13,
     -  0, 2, 3, 4, 5, 6, 6, 7, 8, 8, 9,10,10,11,11,12,12,13,13,13,
     -  0, 2, 3, 4, 5, 6, 6, 7, 8, 9, 9,10,10,11,12,12,13,13,13,14/
c     Maximum critical values for correlation test:
      data nmxcrt/60*0,
     -  0, 0, 0, 0, 9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     -  0, 0, 0, 9,10,10,11,11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     -  0, 0, 0, 9,10,11,12,12,13,13,13,13, 0, 0, 0, 0, 0, 0, 0, 0,
     -  0, 0, 0, 0,11,12,12,12,14,14,14,14,15,15,15, 0, 0, 0, 0, 0,
     -  0, 0, 0, 0,11,12,13,13,14,15,15,16,16,16,16,17,17,17,17,17,
     -  0, 0, 0, 0, 0,13,14,14,15,16,16,16,17,17,18,18,18,18,18,18,
     -  0, 0, 0, 0, 0,13,14,15,16,16,17,17,18,18,18,19,19,19,20,20,
     -  0, 0, 0, 0, 0,13,14,15,16,17,17,18,19,19,19,20,20,20,21,21,
     -  0, 0, 0, 0, 0,13,14,16,16,17,18,19,19,20,20,21,21,21,22,22,
     -  0, 0, 0, 0, 0, 0,15,16,17,18,19,19,20,20,21,21,22,22,23,23,
     -  0, 0, 0, 0, 0, 0,15,16,17,18,19,20,20,21,22,23,23,23,23,24,
     -  0, 0, 0, 0, 0, 0,15,16,18,18,19,20,21,22,22,23,23,24,24,25,
     -  0, 0, 0, 0, 0, 0, 0,17,18,19,29,21,21,22,23,23,24,25,25,25,
     -  0, 0, 0, 0, 0, 0, 0,17,18,19,20,21,22,23,23,24,25,25,26,26,
     -  0, 0, 0, 0, 0, 0, 0,17,18,19,20,21,22,23,24,25,25,26,26,27,
     -  0, 0, 0, 0, 0, 0, 0,17,18,29,21,22,23,23,24,25,26,26,27,27,
     -  0, 0, 0, 0, 0, 0, 0,17,18,20,21,22,23,25,25,26,26,27,27,28/
c     print *,'BATCHMEAN START iout=',iout,' NPTS=',npts
      if (nopr .eq. 0) write (iout,1001) label(1:llabel),npts,nskip
      sd2i=0.0
      ci=999.0
      n=npts-nskip
      if (n .lt. 2) return
c     Calculate avg, sd on the full data set
      datsum=0.d0
      datsum2=0.d0
      do i=1,n
        datsum=datsum+u(i-nskip)
        datsum2=datsum2+u(i-nskip)**2
      end do
      avg=datsum/dfloat(n)
      sd=dsqrt(datsum2/dfloat(n)-(datsum/dfloat(n))**2)
c     print *,'AVG,SD=',avg,sd,' N=',n
c     Aggregate, if needed to max MAXSORT blocks
      if (n .gt. MAXSORT) then
        lb=n/MAXSORT
        nb=MAXSORT
        do ib=1,nb
          datsum=0.d0
          do i=(ib-1)*lb+1,ib*lb
            datsum=datsum+u(i-nskip)
          end do
          dat(ib)=datsum/dfloat(lb)
        end do
        n=nb
      else
        do i=nskip+1,npts
          dat(i-nskip)=u(i)
        end do
      end if
      lenblk=1
      call blankout(decide,1,12)
      do while (n .ge. 4 .and. decide .ne. low)
        call indexit(ixdat,1,n,0)
        call trnsfr(datsort,dat,n)
        call mrgsrt(iout,ixdat,datsort,n,ifst,ilst,itemp,temp,n)
        if (mod(n,2) .eq. 1) then
          datmed=datsort(n/2+1)
        else
          datmed=(datsort(n/2)+datsort(n/2+1))/2.0
        end if
        Z=0.0
        nup=0
        ndown=0
        nrun=1
        rp=dat(1)-datmed
        do i=2,n
          if (rp*(dat(i)-datmed) .le. 0.0) nrun=nrun+1
          if (rp .gt. 0.0) nup=nup+1
          if (rp .le. 0.0) ndown=ndown+1
          rp=dat(i)-datmed
        end do
        if (mod(n,2) .eq. 1) nrun=nrun-1
        if (rp .gt. 0.0) nup=nup+1
        if (rp .le. 0.0) ndown=ndown+1
        if (ndown .eq. 0 .or. nup .eq. 0) then
          decide=low
        else if (max0(nup,ndown) .le. 20) then
          if (nmxcrt(ndown,nup) .eq. 0 .or.
     -        nmncrt(ndown,nup) .eq. 0) then
            decide=low
          else
c           Use the table
            if (nrun .lt. nmncrt(ndown,nup) .or.
     -           nrun .gt. nmxcrt(ndown,nup)) then
              decide=corr
            else
              decide=uncorr
            end if
          end if
        else
c         Use formula assuming normal distribution
          r=float(2*nup*ndown)/float(nup+ndown)+1
          s2=float(2*nup*ndown*(2*nup*ndown-nup-ndown))/
     -      float((nup+ndown)**2*(nup+ndown))
          z=(float(nrun)-r)/sqrt(s2)
          if (z .gt. 1.96) then
            decide=corr
          else
            decide=uncorr
          end if
        end if
        if (decide .eq. uncorr) then
          datsum=0.d0
          datsum2=0.d0
          do i=1,n
            datsum=datsum+dat(i)
            datsum2=datsum2+dat(i)**2
          end do
          avg=datsum/dfloat(n)
          sd=dsqrt(datsum2/dfloat(n)-(datsum/dfloat(n))**2)
          ci=sd/sqrt(float(n-1))    
        end if
        if (nopr .eq. 0) write (iout,1002)
     -    label(1:llabel),lenblk,nup,ndown,nrun,decide,ci,z
        if (ci .ne. 999.0) return
        do i=1,n/2
          dat(i)=(dat(2*i)+dat(2*i-1))/2.0
        end do
        n=n/2
        lenblk=2*lenblk
      end do
      return
1001  format(/,1x,a8,' Number of data points=',i6,
     -  ' Number of data points skipped=',i4)
1002    format(1x,a,' block size=',i9,' nup=',i3,' ndown=',i3,
     -    ' nrun=',i3,2x,a12,' CI(1s)=',f10.4,' Z=',f10.4)
      end
