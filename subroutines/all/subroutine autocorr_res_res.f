      subroutine autocorr_res_res(nbres,nframe,iframeunit,
     -  framefac,itrack,nhbdist,nhbpers,maxlenon,maxlenoff,
     -  itrackf,itrackl,ifres,iresno,isegno,resnames,nrescol,
     -  it1,it2,iuselaston,iout,mxbonds,mxframes,mxrsd,mxrec)
      dimension itrack(mxframes),
     -  nhbdist(mxbonds),nhbpers(mxbonds),maxlenon(mxbonds),
     -  maxlenoff(mxbonds),itrackf(mxbonds),itrackl(mxbonds),
     -  ifres(mxrsd),iresno(mxrec),isegno(mxrec),it1(mxbonds),
     -  it2(mxbonds)
      character*8 resnames(mxrsd)
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (6*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
     -  ihbtores(MAXBONDS),iusepair(MAXBONDS),nhb_atot(MAXBONDS),
     -  nhb_rtot(MAXBONDS),a1(MAXBONDS),a2(MAXBONDS),fill(IFILL2)
      common /bondpairs/ ihbpair(2,MAXBONDS),ihb_pair_res(3,MAXBONDS)
      parameter (MAXFRAMES=50000,MAXCOPY=600,MAXITEMS=2*MAXCOPY-2)
c     print *,'AUTOCORR_RES_RES NRESCOL=',nrescol,
c    -  ' MXRSD,MXREC=',mxrsd,mxrec
      nreusemax=nframe
      call auc_params(iauctype,lastframeinp,loffmin,nreusemax,
     -  nseg_scr,iaucw,nframe,iframeunit,framefac,'Residue-residue',15,
     -  iuselaston,iout,mxframes)
      call trackstat(nbres,nhbdist,nhbpers,maxlenon,maxlenoff,itrackf,
     -  itrackl,it1,it2)
      nauc_extra=0
      do irr=1,nbres
c       Get the irr-th track
        ir1=ihb_pair_res(1,irr)
        ir2=ihb_pair_res(2,irr)
        perc=100.0*float(ihb_pair_res(3,irr))/float(nframe)
        ia1=ifres(ir1)
        ia2=ifres(ir2)
c       write (iout,2000) irr,
c    -    line(index(ia1))(irescol1:irescol2),iresno(ia1),isegno(ia1),
c    -    line(index(ia2))(irescol1:irescol2),iresno(ia2),isegno(ia2),perc
        write (iout,2000) irr,
     -    resnames(ir1)(1:nrescol),iresno(ia1),isegno(ia1),
     -    resnames(ir2)(1:nrescol),iresno(ia2),isegno(ia2),perc
        call persistence(nhbdist(irr),nhbpers(irr),itrackf(irr),
     -    itrackl(irr),maxlenon(irr),maxlenoff(irr),iframeunit,
     -    framefac,nframe,iuselaston,iout)
        call getbondtrack(irr,itrack,ifirstframe,lastframe,30,nframe)
        call autocorr(irr,irr,itrack,ifirstframe,lastframe,iframeunit,
     -    framefac,iauctype,lastframeinp,nreusemax,percon,loffmin,
     -    nseg_scr,nauc_extra,nframe,iout,mxframes)
      end do
      return
2000  format(' RR#',i4,1x,a,i5,' C/S',i2,' - ',a,i5,' C/S',i2,
     -  ' perc=',f5.1)
      end
