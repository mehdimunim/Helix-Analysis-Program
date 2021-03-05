      subroutine histogram(r,n,rmaxi,rbin,label,llabel,ihistogram,iout,
     -  mxh)
      dimension r(n)
      dimension ihistogram(mxh)
      character*(*) label
      if (rmaxi .eq. 0.0) then
        rmin=1.e32
        rmax=-rmin
        do i=1,n
          if (r(i) .lt. rmin) rmin=r(i)
          if (r(i) .gt. rmax) rmax=r(i)
        end do
      else
        rmin=0.0
        rmax=rmaxi
      end if
      ng=(rmax-rmin)/rbin+1.0
      if (ng .lt. 2) then
        write (6,1001) rmin,rmax,rbin
        return
      end if
      call zeroiti(ihistogram,0,ng)
      do i=1,n
        ix=ng*(r(i)-rmin)/(rmax-rmin)+1 
        if (ix .lt. 1) ix=1
        if (ix .gt. ng) ix=ng
        ihistogram(ix)=ihistogram(ix)+1
      end do
      write (iout,1008) n,label(1:llabel),rbin,rmin,rmax,
     -  (ihistogram(i),i=1,ng)
c     imax=0
c     do i=1,ng
c       if (ihistogram(i) .gt. imax) imax=ihistogram(i)
c     end do
c     if (imax .gt. ng) then
c       do i=1,ng
c         ihistogram(i)=10.0*float(ihistogram(i))/float(imax)
c       end do
c     end if
      return
1001  format(' Range: [',e12.5,' - ',e12.5,'] and bin size ',f8.3,
     -  ' does not give a histogram') 
1008  format(' Histogram of ',i5,1x,a,/,' Binsize=',f8.3,' in the ',
     -  f8.3,' - ',f8.3,' range:',/,(10i6))
      end
