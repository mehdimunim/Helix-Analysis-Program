      subroutine trackstat(nbonds,nhbdist,nhbpers,maxlenon,
     -  maxlenoff,itf,itl,it,itprev)
      dimension nhbdist(nbonds),nhbpers(nbonds),maxlenon(nbonds),
     -  maxlenoff(nbonds),itf(nbonds),itl(nbonds),it(nbonds),
     -  itprev(nbonds)
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (6*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
     -  ihbtores(MAXBONDS),nusepair(MAXBONDS),nhb_atot(MAXBONDS),
     -  nhb_rtot(MAXBONDS),a1(MAXBONDS),a2(MAXBONDS),fill(IFILL2)
      parameter (MAXFRAMES=50000,MAXCOPY=600,MAXITEMS=2*MAXCOPY-2)
      common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
     -  increment,increment2,ires(MAXITEMS,MAXFRAMES),
     -  scres(2,MAXFRAMES),x0(MAXCOPY),y0(MAXCOPY),nxselres,
     -  ixselres(MAXCOPY)
      dimension isegstart(MAXBONDS)
c     print *,'TRACKSTAT nbonds=',nbonds,' nframe=',nframe
      call zeroiti(nhbdist,0,nbonds)
      call zeroiti(nhbpers,0,nbonds)
      call zeroiti(maxlenon,0,nbonds)
      call zeroiti(maxlenoff,0,nbonds)
      call zeroiti(isegstart,0,nbonds)
      call zeroiti(itf,0,nbonds)
      call zeroiti(itl,0,nbonds)
      do ifr=1,nframe
        call readbitc(ires(1,ifr),it,nbonds,30,MAXITEMS)
        do i=1,nbonds
          if (itf(i) .gt. 0) then
            if (it(i) .ne. itprev(i)) then
              len=ifr-isegstart(i)
               if (itprev(i) .eq. 1) then
                 if (maxlenon(i) .lt. len) maxlenon(i)=len
               else
                 if (maxlenoff(i) .lt. len) maxlenoff(i)=len
               end if
               isegstart(i)=ifr
            end if
          end if
          if (it(i) .eq. 1) then
            nhbdist(i)=nhbdist(i)+1
            if (itf(i) .eq. 0) then
              itf(i)=ifr
              isegstart(i)=ifr
            end if
            itl(i)=ifr
            if (ifr .gt. 1 .and. itprev(i) .eq. 0)
     -        nhbpers(i)=nhbpers(i)+1
          end if
        end do
        call trnsfi(itprev,it,nbonds)
      end do
      return
      end
