      subroutine rmsf_av(cdp,cdp2,atw,n,index,ixres,rmsfsum,rmsfav,
     -  nframe,nslt,nresslt)
      real*8 cdp(3,n),cdp2(3,n),rmsfsum(nresslt)
      dimension atw(nslt),index(n),ixres(nslt),rmsfav(nresslt)
      real*8 atwrsum
      call zeroitd(rmsfsum,nresslt)
      ia=1
      ir=ixres(index(ia))
      do while (ia .le. n)
        rmsfs=0.0
        atwrsum=0.d0
        do while (ir .eq. ixres(index(ia)))
c          write (77,9671) ia,ir,index(ia),atw(index(ia)),
c     -      (cdp(k,ia),k=1,3),(cdp2(k,ia),k=1,3)
c9671      format(i4,' ir=',i3,' ix=',i4,' atw=',f8.3,' cdp=',3e12.5,
c     -      ' cdp2=',3e12.5)
          atwrsum=atwrsum+atw(index(ia))
          do k=1,3
            rmsfsum(ir)=rmsfsum(ir)+
     -        atw(index(ia))*(cdp2(k,ia)/nframe-(cdp(k,ia)/nframe)**2)
          end do
          ia=ia+1
          if (ia .gt. n) go to 100
        end do
100     continue
        rmsfav(ir)=dsqrt(rmsfsum(ir)/atwrsum)
c       write (77,*) 'ir=',ir,' RMSFAV=',rmsfav(ir)
        if (ia .le. n) ir=ixres(index(ia))
      end do
      return
      end
