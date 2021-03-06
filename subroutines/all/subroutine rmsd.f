      subroutine rmsd(co,c,n,nfinalov,nfinalrmsd,atw,atwsum,atw1,c1,
     -  c2,indxov1,indxov2,noopt2d,isymm,indxrmsd1,indxrmsd2,rot,com1,
     -  com2,etot,etot2,iout,devmax,devmaxnoopt,maxat)
      dimension co(3,n),c(3,n),atw(n),atw1(n),c1(3,n),c2(3,n),
     -  indxov1(n),indxov2(n),indxrmsd1(n),indxrmsd2(n),rot(3,3),
     -  com1(3),com2(3)
      real*8 atwsum
c     Calculate the RMSD of the two structures after obtaining the best fit
c     using Kabsch formula
      character*2 iatnm2
      common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
     -  mmatno(64),iatnm2(99)
      parameter (MAXFRAMES=50000,MAXCOPY=600)
      common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
     -  increment,increment2,res(2,MAXFRAMES,MAXCOPY),
     -  x0(MAXCOPY),y0(MAXCOPY),nxselres,ixselres(MAXCOPY)
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (4*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
     -  ihbpair(2,MAXBONDS),a1(MAXBONDS),a2(MAXBONDS),fill(IFILL2)
c     print *,'Start RMSD nframeref,nframe=',nframeref,nframe,
c    -  ' nfinalov,nfinalrmsd=',nfinalov,nfinalrmsd
c     Calculate rMSD w/o overlay
c      write (iout,8811) n,nfinalrmsd,nfinalov,iedit,limresrange
c8811  format(' n,nfinalrmsd,nfinalov,iedit,limresrange=',3i8,2i3)
c      write (iout,8812) 'indxov1',(indxov1(i),i=1,nfinalov)
c      write (iout,8812) 'indxrmsd1',(indxov2(i),i=1,nfinalrmsd)
c8812  format(1x,a,':',/,(20i4))
c      write (iout,8813) 'C',((c(k,indxov1(i)),k=1,3),i=1,nfinalov)
c      write (iout,8813) 'CO',((co(k,indxov2(i)),k=1,3),i=1,nfinalov)
c8813  format(1x,a,':',/,(5f10.5))
c      write (iout,8814) (atw(indxov1(i)),i=1,nfinalrmsd)
c8814  format(' ATW:',/,(5f10.5))
      sdnoopt=sqrt(sdsum(nfinalrmsd,co,c,atw,indxrmsd1,indxrmsd2,
     -  devmaxnoopt,maxat))
      if (noopt2d .eq. 0) then
        call bestoverlay(nfinalov,indxov1,indxov2,co,c,atw,atwsum,c1,c2,
     -    atw1,rot,com1,com2,00,0.001,iout,maxat)
c       Shift the full set of coordinate to com1,com2, and rotate by rot
        call shiftmol(co,n,com1,c1,-1.0)
        call shiftmol(c,n,com2,c2,-1.0)
        call rotate_c(c2,n,rot,c2,'RMSD',4)
        sd=sqrt(sdsum(nfinalrmsd,c1,c2,atw,indxrmsd1,indxrmsd2,devmax,
     -    maxat))
      else
        devmax=devmaxnoopt
      end if
c      do i=1,45
c        write (78,8674)i,(co(k,i),k=1,3),(c1(k,i),k=1,3),(c2(k,i),k=1,3)
c8674    format(i3,' co=',3f10.5,' c1=',3f10.5,' c2=',3f10.5)
c      end do
      call trajlimtest(nframe,MAXFRAMES)
      if (etot .eq. 0.0) then
        if (noopt2d .eq. 0)
     -    write (iout,1001) 'Best overlap:',sd,devmax
        write (iout,1001) 'No   overlap:',sdnoopt,devmaxnoopt
      else
        if (noopt2d .eq. 0) then
          write (iout,1001) 'Best overlap:',sd,devmax,' ','c',etot
          write (iout,1001) 'No   overlap:',sdnoopt,devmaxnoopt,' ',
     -      'r',etot2
        else
          write (iout,1001) 'No   overlap:',sdnoopt,devmaxnoopt,' ',
     -      'c',etot
        end if
        res(1,nframe,9)=etot
      end if
      if (nframeref .le. 1) then
        res(1,nframe,7)=sd
        res(2,nframe,7)=devmax
        res(1,nframe,8)=sdnoopt
        res(2,nframe,8)=devmaxnoopt
      end if
      if (nframeref .gt. 0) then
        if (noopt2d .eq. 0) then
          rmsd2d(nframe,nframeref)=sd
          if (isymm .eq. 1) rmsd2d(nframeref,nframe)=sd
        else
          rmsd2d(nframe,nframeref)=sdnoopt
          if (isymm .eq. 1) rmsd2d(nframeref,nframe)=sdnoopt
        end if
      end if
      return
1001  format(1x,a,' RMSD=',f8.3,' Max dev=',f10.3,' A',
     -  a,'E',a1,'=',e14.7)
      end
