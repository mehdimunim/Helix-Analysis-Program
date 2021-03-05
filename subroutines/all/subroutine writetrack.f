      subroutine writetrack(iout_track,iout,nbits,ntracks,nres2d,
     -  trackfile,ltrackfile,ianc_anc)
      dimension ianc_anc(ntracks)
      character*(*) trackfile
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (6*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
     -  ihbtores(MAXBONDS),nusepair(MAXBONDS),nhb_atot(MAXBONDS),
     -  nhb_rtot(MAXBONDS),a1(MAXBONDS),a2(MAXBONDS),fill(IFILL2)
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
      write (iout_track,1000) trajnam(1:ltrajnam)
      write (iout_track,1001)nframe,ntracks,nres2d,firsttraj,ilasttraj,
     -  incrementtraj
      do itr=1,ntracks
        call zeroiti(itrack,0,nframe)
        ii=(itr-1)/nbits+1
        ib=mod(itr-1,nbits)
        do ifr=1,nframe
          if (btest(ires(ii,ifr),ib)) itrack(ifr)=1
        end do
        write (iout_track,1006) itr,(ihbpair(k,itr),k=1,2),
     -    ianc_anc(itr),(itrack(ifr),ifr=1,nframe)
      end do
      write (iout,1002) trackfile(1:ltrackfile),nframe,ntracks
      close (iout_track)
      return
1000  format(a)
1001  format(6i7)
1002  format(' Bond tracks were written to file ',a,/,
     -  ' Number of frames=',i7,/,' Number of tracks=',i5)
1006  format(' TRACK#',i5,' ia1=',i7,' ia2=',i7,' iAA=',i1,/,(100i1))
      end
