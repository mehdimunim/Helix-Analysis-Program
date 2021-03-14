      subroutine arminmax2(ar,nfst,nlst,nplot,armin1,armax1,armin2,
     -  armax2,iyinc,ndim)
      dimension ar(ndim,nlst)
c     if iyinc=1, use the 2nd column of ar (valid only for nplot=1)
c     print *,'ARMINMAX2 nfst,nlst=',nfst,nlst
      if (nplot .eq. 2 .and. iyinc .gt. 0) then
        write (6,2000) nplot,iyinc
        iyinc=0
      end if
      armin1=1.e+30
      armax1=-armin1
      do i=nfst,nlst
        if (ar(iyinc+1,i) .lt. armin1) armin1=ar(iyinc+1,i)
        if (ar(iyinc+1,i) .gt. armax1) armax1=ar(iyinc+1,i)
      end do
      if (nplot .eq. 2) then
        armin2=1.e+30
        armax2=-armin2
        do i=nfst,nlst
          if (ar(2,i) .lt. armin2) armin2=ar(2,i)
          if (ar(2,i) .gt. armax2) armax2=ar(2,i)
        end do
      end if
      return
2000  format(' PROGRAM ERROR: illegal nplot,iyinc (',i2,',',i2,')')
      end
