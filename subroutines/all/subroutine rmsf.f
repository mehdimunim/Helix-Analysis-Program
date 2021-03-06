      subroutine rmsf(cref,c,atw,nslt,ixres,index,bfacavg,nresslt,
     -  maxat)
      dimension cref(3,maxat),c(3,maxat),atw(maxat),ixres(maxat),
     -  index(maxat)
      real*8 bfacavg(nresslt)
c     print *,'RMSF nslt=',nslt
c      write (78,9878) (i,index(i),index(i),atw(index(i)),
c     -  (cref(k,index(i)),k=1,3),(c(k,index(i)),k=1,3),i=1,nslt)
c9878  format(i5,' ix1,2=',2i4,' aw=',f8.3,' c1=',3f10.5,' c2=',3f10.5)
c     write (40,8789) (index(ia),ia=1,nslt)
c8789 format(20i4)
      ia=1
      ir=ixres(index(ia))
      do while (ia .le. nslt)
        rmsfs=0.0
        atwrsum=0.0
        do while (ir .eq. ixres(index(ia)))
          iaa=index(ia)
          dev2=0
          do k=1,3
            dev2=dev2+(c(k,iaa)-cref(k,iaa))**2
          end do
          rmsfs=rmsfs+atw(iaa)*dev2
          atwrsum=atwrsum+atw(iaa)
c         write (76,9782) ia,iaa,ir,atw(iaa),(cref(k,iaa),k=1,3)
c9782     format(' ia,index(ia),ir=',3i5,' aw=',f8.3,' cref=',3f10.5)
          ia=ia+1
          if (ia .gt. nslt) go to 100
          if (index(ia) .eq. 0) print *,'INDEX(',ia,')=0 !!'
        end do
100     continue
c       bfacavg(ir)=bfacavg(ir)+sqrt(rmsfs/atwrsum)
        bfacavg(ir)=bfacavg(ir)+rmsfs/atwrsum
        if (ia .le. nslt) ir=ixres(index(ia))
      end do
      return
      end
