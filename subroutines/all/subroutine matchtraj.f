      subroutine matchtraj(rmsdsim,iw0)
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (4*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
     -  ihbpair(2,MAXBONDS),a1(MAXBONDS),a2(MAXBONDS),fill(IFILL2)
      parameter (MAXFRAMES=50000,MAXCOPY=600)
      common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
     -  increment,increment2,res(2,MAXFRAMES,MAXCOPY),
     -  x0(MAXCOPY),y0(MAXCOPY),nxselres,ixselres(MAXCOPY)
      dimension nsame1(MAX2D),nsame2(MAX2D)
      character*6 lab
      nframe1=nframeref
      nframe2=nframe
      call zeroiti(nsame1,0,nframe1)
      call zeroiti(nsame2,0,nframe2)
      rmsdmin=100000.0
c     print *,'MATCHTRAJ iw0,nframe,nframeref=',iw0,nframe,nframeref
      do if=1,nframe1
c       Find closest to frame if
        rmsdmini=100000.0
        do jf=1,nframe2
          if (rmsdmini .gt. rmsd2d(if,jf)) then
            rmsdmini=rmsd2d(if,jf)
            jmin=jf
          end if
          if (rmsdsim .gt. rmsd2d(if,jf)) then
            nsame1(if)=nsame1(if)+1
          end if
        end do
        a1(if)=jmin
        if (rmsdmin .gt. rmsdmini) then
          rmsdmin=rmsdmini
          igmin=if
          jgmin=jmin
        end if
      end do
      do jf=1,nframe2
c       Find closest to frame if
        rmsdminj=100000.0
        do if=1,nframe1
          if (rmsdminj .gt. rmsd2d(if,jf)) then
            rmsdminj=rmsd2d(if,jf)
            imin=if
          end if
          if (rmsdsim .gt. rmsd2d(if,jf)) then
            nsame2(jf)=nsame2(jf)+1
          end if
        end do
        a2(jf)=imin
      end do
      write (iw0,1004) 2,1
      do if=1,nframe1
        jmin=a1(if)
        minrev1=a2(jmin)
        lab='      '
        if (minrev1 .eq. if) lab='mutual'
        if (rmsdsim .eq. 0) then
          write (iw0,1000) 1,if,2,jmin,rmsd2d(if,jmin),lab
        else
          write (iw0,1003) 1,if,rmsdsim,nsame1(if),jmin,
     -      rmsd2d(if,jmin),lab
        end if
      end do
      write (iw0,1004) 1,2
      do jf=1,nframe2
        imin=a2(jf)
        minrev2=a1(imin)
        lab='      '
        if (minrev2 .eq. jf) lab='mutual'
        if (rmsdsim .eq. 0) then
          write (iw0,1000) 2,jf,1,imin,rmsd2d(imin,jf),lab
        else
          write (iw0,1003) 2,jf,rmsdsim,nsame2(jf),imin,
     -      rmsd2d(imin,jf),lab
        end if
      end do
      write (iw0,1001) igmin,jgmin,rmsdmin
      return
1000  format(' Traj',i2,' frame',i5,' closest from traj',i2,' is frame',
     -  i5,' RMSD=',f10.1,1x,a)
1001  format(' Best match is between traj 1, frame',i5,' traj 2, frame',
     -  i5,' RMSD=',f10.1)
1003  format(' Tr',i2,' Fr',i5,' N(<',f4.1,')=',i2,' closest:',i5,
     -  ' RMSD=',f5.1,1x,a)
1004  format(' Matching Trajectory',i2,' frames to Trajectory',i2,
     -  ' frames')
      end
