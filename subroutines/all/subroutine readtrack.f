      subroutine readtrack(iout_track,iout,nbits,nbfound,nbresfound,
     -  nres2d,ixres,trackfile,ltrackfile,ianc_anc,ifail,maxrec)
      dimension ixres(maxrec),ianc_anc(nbfound)
      character*(*) trackfile
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (6*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
     -  ihbtores(MAXBONDS),nusepair(MAXBONDS),nhb_atot(MAXBONDS),
     -  nhb_rtot(MAXBONDS),itempres(MAXBONDS),a2(MAXBONDS),fill(IFILL2)
      common /bondpairs/ ihbpair(2,MAXBONDS),ihb_pair_res(3,MAXBONDS)
      parameter (MAXFRAMES=50000,MAXCOPY=600,MAXITEMS=2*MAXCOPY-2)
      common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
     -  increment,increment2,ires(MAXITEMS,MAXFRAMES),
     -  scres(2,MAXFRAMES),x0(MAXCOPY),y0(MAXCOPY),nxselres,
     -  ixselres(MAXCOPY)
      character*200 trajnam,trajnam2
      common /trajname/ trajnam,trajnam2,ltrajnam,ltrajnam2,ifirsttraj,
     -  ifirsttraj2,ilasttraj,ilasttraj2,incrementtraj,incrementtraj2
      dimension itrack(MAXFRAMES)
      itr=0
      call blankout(trajnam,1,200)
      read (iout_track,1000,end=888,err=888) trajnam
      call lastchar(trajnam,ltrajnam,200)
      nres2d=0
      read (iout_track,1001,end=100,err=100) nframe,ntracks,nres2d,
     -  ifirsttraj,ilasttraj,incrementtraj
100   if (nres2d .eq. 0) go to 888
      nframetot=nframe
      nbfound=ntracks
      do itr=1,ntracks
        call zeroiti(itrack,0,nframe)
        read (iout_track,1004,end=888,err=888) itrr,
     -    (ihbpair(k,itr),k=1,2),ianc_anc(itr),
     -    (itrack(ifr),ifr=1,nframe)
        ii=(itr-1)/nbits+1
        ib=mod(itr-1,nbits)
        do ifr=1,nframe
          if (itrack(ifr) .eq. 1)
     -      ires(ii,ifr)=ibset(ires(ii,ifr),ib)
        end do
        ir1=ixres(ihbpair(1,itr))
        ir2=ixres(ihbpair(2,itr))
        ihb=1
        do while (ihb .le. nbresfound .and.
     -    (ir1 .ne. ihb_pair_res(1,ihb) .or.
     -     ir2 .ne. ihb_pair_res(2,ihb)))
          ihb=ihb+1
        end do
        if (ihb .gt. nbresfound) then
          nbresfound=nbresfound+1
          ihb_pair_res(1,nbresfound)=ir1
          ihb_pair_res(2,nbresfound)=ir2
          ihb_pair_res(3,nbresfound)=1
          ihb_pair_res(3,nbresfound+1)=0
          itempres(nbresfound)=1
        else if (itempres(ihb) .eq. 0) then
          ihb_pair_res(3,ihb)=ihb_pair_res(3,ihb)+1
          itempres(ihb)=1
        end if
      end do
      write (iout,1002) trackfile(1:ltrackfile)
      write (iout,1003) nframe,ntracks,nbresfound
      write (6,1003) nframe,ntracks,nbresfound
      call write_traj_lim(iout,' ',0,1,incr_tr,0)
      close (iout_track)
      ifail=0
      return
888   if (itr .eq. 0) write (6,1005)
      if (itr .gt. 0) write (6,1006) itr
      return
1000  format(a)
1001  format(6i7)
1002  format(' Bond tracks were read from file ',a)
1003  format(' Number of frames=',i7,/,' Number of tracks=',i5,/,
     -  ' Number of residue-residue pairs=',i4)
1004  format(' TRACK#',i5,' ia1=',i7,' ia2=',i7,' iAA=',i1,/,(100i1))
1005  format(' ERROR reading the track file header')
1006  format(' ERROR reading track # ',i5)
      end
