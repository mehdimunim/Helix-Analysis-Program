      subroutine read_2drmsd(iw0,system,lsystem,trajname,ltrajname,
     -  trajname2,ltrajname2,ifirst,ilast,incr,ifirst2,ilast2,
     -  incr2,isym,rmsdmn,rmsdmx,nframeread,nframex,nframey,
     -  ietotsaved,maxecho,limresrange,ierr)
      character*(*) trajname,trajname2,system
      parameter (MAXFRAMES=50000,MAXCOPY=600)
      parameter (MAXCOPY1=MAXCOPY-1)
      common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
     -  increment,increment2,res(2,MAXFRAMES,MAXCOPY1),
     -  xyplot(2,MAXFRAMES),x0(MAXCOPY),y0(MAXCOPY),nxselres,
     -  ixselres(MAXCOPY)
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (4*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
     -  ihbpair(2,MAXBONDS),a1(MAXBONDS),a2(MAXBONDS),fill(IFILL2)
      character*16 errlab
      character*80 line
      ierr=1
      limresrange=1
      nline=0
      rmsdmn=1000000.0
      ltrajname2=0
      if (isym .eq. 1) rmsdmn=0.0
      rmsdmx=0.0
      nframex=0
      nframey=0
      call blankout(errlab,1,16)
      do while (.true.)
        call blankout(line,1,80)
        read (iw0,1000,end=999) line
        nline=nline+1
        if (line(1:9) .eq. ' Config #') then
c         RMSD data found
          errlab='structure number'
cConfig #     1 (inp #:          1) Ref #     1 (inp #:         1)
          read (line(11:15),*,err=888) ix
          read (line(43:47),*,err=888) jx
c         read (line(52:57),*,err=888) jx
          if (ix .gt. nframex) nframex=ix
          if (jx .gt. nframey) nframey=jx
          read (iw0,1000,end=999) line
          errlab='RMSD'
          read (line(21:28),*,err=888) rmsd2d(ix,jx)
          if (rmsd2d(ix,jx) .gt. rmsdmx) rmsdmx=rmsd2d(ix,jx)
          if (rmsd2d(ix,jx) .lt. rmsdmn) rmsdmn=rmsd2d(ix,jx)
          if (isym .eq. 1) rmsd2d(jx,ix)=rmsd2d(ix,jx)
          if (line(51:53) .eq. 'Ec=') then
            errlab='energy'
            read (line(64:77),*,err=888) res(1,ix,11)
            ietotsaved=1
          end if
          nline=nline+1
          lc=2
          do while (lc .gt. 1)
            call blankout(line,1,80)
            read (iw0,1000,end=999) line
            call lastchar(line,lc,80)
            nline=nline+1
          end do
        else if (line(1:33) .eq. ' Second trajectory file analyzed:')
     -           then
          call lastchar(line,icl,80)
          ltrajname2=icl-33
          trajname2(1:ltrajname2)=line(34:33+ltrajname2)
          nline=nline+1
          print *,'Second trajectory: ',trajname2(1:ltrajname2)
        else if (line(1:26) .eq. ' Trajectory file analyzed:') then
          call lastchar(line,icl,80)
          ltrajname=icl-26
          trajname(1:ltrajname)=line(27:26+ltrajname)
          nline=nline+1
          print *,'Trajectory: ',trajname(1:ltrajname)
        else if (line(1:47) .eq.
     -    ' Trajectory first and last frame and increment:') then
          read (line(48:68),1001,end=777,err=777) ifirst,ilast,incr
          call save_traj_lim(ifirst,ilast,incr,1)
         else if (line(1:54) .eq.
     -    ' Second trajectory first and last frame and increment:') then
          read (line(55:75),1001,end=777,err=777) ifirst2,ilast2,incr2
          call save_traj_lim(ifirst2,ilast2,incr2,2)
        else if (line(1:27) .eq. ' Number of frames analyzed=') then
          call lastchar(line,icl,80)
          errlab='number of frames'
          read (line(28:icl),*,err=888) nframeread
          nline=nline+1
          print *,'Number of structures analyzed=',nframeread
          go to 100
        else if (line(1:8) .eq. ' System:') then
          call blankout(line,1,80)
          read (iw0,1000,end=999) line
          call lastchar(line,lsystem,80)
          system(1:lsystem)=line(1:lsystem)
          nline=nline+2
          print *,'System:',system(1:lsystem)
        else if (line(1:28) .eq. ' All atoms are used for RMSD') then
          limresrange=0
        end if
      end do
100   if (ltrajname2 .gt. 0) then
c       Cross RMSD map
        do i=1,min0(maxecho,nframex)
          write (iw0,2000) i,(rmsd2d(i,j),j=1,nframey)
        end do
      else
        do i=2,min0(maxecho,nframeread)
          write (iw0,2000) i,(rmsd2d(j,i),j=1,i-1)
        end do
        nframex=max0(nframex,nframey)
        nframey=max0(nframex,nframey)
      end if
      write (6,2001) nframex,nframey
c     write (78,2001) nframex,nframey
c     do i=1,10
c       write (78,*) 'RMSFCLUSTER i=',i,'RMSD:'
c       write (78,8712) (rmsd2d(i,j),j=1,nframex)
c8712   format(10f8.4)
c     end do
      ierr=0
      return
777   print *,'ERROR: invalid frame limits:'
      print *,line
      return
888   print *,'ERROR: invalid input for ',errlab,':'
      print *,line
      return
999   print *,'ERROR: 2D RMSD file ended prematurely at line #',nline
      return
1000  format(a)
1001  format(3i7)
2000  format(' RMSD(',i4,'):',10f6.1,/,(12x,10f6.1))
2001  format(' Matrix of ',i6,' rows and ',i6,' columns read')
      end
