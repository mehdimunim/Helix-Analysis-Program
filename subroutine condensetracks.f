      subroutine condensetracks(nbfoundorig,nbres,ibframe,ibframe_rr,
     -  iout,mxbonds)
      dimension ibframe(mxbonds),ibframe_rr(mxbonds)
      parameter (MAXFRAMES=50000,MAXCOPY=600,MAXITEMS=2*MAXCOPY-2)
      common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
     -  increment,increment2,ires(MAXITEMS,MAXFRAMES),
     -  scres(2,MAXFRAMES),x0(MAXCOPY),y0(MAXCOPY),nxselres,
     -  ixselres(MAXCOPY)
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (6*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
     -  ihbtores(MAXBONDS),nusepair(MAXBONDS),nhb_atot(MAXBONDS),
     -  nhb_rtot(MAXBONDS),ibond(MAXBONDS),a1(MAXBONDS),fill(IFILL2)
c     Condense the bond-bits to res-res bond bits in ires(*,nframe)
      write (iout,2000) nbfoundorig,nbres
      do ifr=1,nframe
        call zeroiti(ibframe_rr,0,nbres)
        call readbitc(ires(1,ifr),ibframe,nbfoundorig,30,MAXITEMS)
        do ib=1,nbfoundorig
          ibb=ihbtores(ib)
          if (ibframe(ib)*ibb .gt. 0) ibframe_rr(ibb)=1
        end do
        call savebitc(ires(1,ifr),ibframe_rr,nbres,30,MAXITEMS)
      end do
      return
2000  format(' Bond tracks of ',i4,' bonds condensed to ',i3,
     - ' res-res tracks')
      end
