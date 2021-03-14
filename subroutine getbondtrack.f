      subroutine getbondtrack(itr,itrack,ifirstframe,lastframe,nbits,
     -  nfrm)
      dimension itrack(nfrm)
      parameter (MAXFRAMES=50000,MAXCOPY=600,MAXITEMS=2*MAXCOPY-2)
      common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
     -  increment,increment2,ires(MAXITEMS,MAXFRAMES),
     -  scres(2,MAXFRAMES),x0(MAXCOPY),y0(MAXCOPY),nxselres,
     -  ixselres(MAXCOPY)
c     Read bits
c     print *,'GETBONDTRACK itr,nbits,nfrm,nfram=',itr,nbits,nfrm,nframe
      call zeroiti(itrack,0,nframe)
      ifirstframe=0
      lastframe=0
      ii=(itr-1)/nbits+1
      ib=mod(itr-1,nbits)
      do ifr=1,nframe
        if (btest(ires(ii,ifr),ib)) itrack(ifr)=1
        if (itrack(ifr) .eq. 1) then
          if (ifirstframe .eq. 0) ifirstframe=ifr
          lastframe=ifr
          if (ifirstframe .eq. 0) ifirstframe=ifr
        end if
      end do
      return
      end
