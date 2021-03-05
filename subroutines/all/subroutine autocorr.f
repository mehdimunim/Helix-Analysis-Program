      subroutine autocorr(itr_ix,itr,itrack,ifirstframe,lastframe,
     -  iframeunit,framefac,iauctypeinp,lastframeinp,nreusemax,percon,
     -  loffmin,nseg_scr,nauc_extra,nframe,iout,mxframes)
      dimension itrack(mxframes)
      parameter (MAXFRAMES=50000)
      parameter (MAXFRAMES2=MAXFRAMES/2)
      dimension auc(MAXFRAMES2),nauc(MAXFRAMES2),lfrunit(4)
      parameter (MAXFRAMES21=MAXFRAMES2+1,MAXCOPY=600)
      common /aucw/ auc_all(MAXFRAMES21,MAXCOPY)
      dimension ran(1),ixran(100)
      character*6 frunit(4)
      data frunit /'frames','ps    ','ns    ','ms    '/,lfrunit/6,2,2,2/
c     Calculate the autocorrelation of bond itr
c     write (6,*) 'AUTOCORR itr_ix,itr,ifirstframe,lastframe=',
c    -  itr_ix,itr,ifirstframe,lastframe
      iauctype=iauctypeinp
      nframe2=max0(nframe,lastframeinp)/2
      call zeroit(auc,nframe2)
      call zeroiti(nauc,0,nframe2)
c      write (0006,9867) itr,nframe,lastframe,loffmin
c9867  format(' ITR=',I6,' NFRAME=',I6,' LASTFR=',I6,' LOFFMIN=',I6)
c     Accumulate auc
      if (nframe-lastframe .gt. loffmin .and. iauctype .gt. 3) then
        iauctype=3
        write (iout,1005) nframe,lastframe,loffmin
      end if
      nframeuse=lastframe-ifirstframe+1
      lastframeuse=0
      if (iauctype .eq. 1) then
        lastframeuse=lastframe
      else if (iauctype .eq. 2) then
        lastframeuse=nframe
      else
        if (iauctype .eq. 3) then
          nframeuse=lastframeinp-ifirstframe+1
          call zeroiti(itrack,nframe,lastframeinp)
        else if (iauctype .eq. 4) then
c         Pad with repeating the track
          nframeuse=lastframeinp-ifirstframe+1
          nadd=lastframeinp-lastframe
c          write (iout,*) 'NADD=',nadd
          incr=lastframe
          lentrackuse=min0(lastframe-ifirstframe+1,nreusemax)
          do while (nadd .gt. 0)
            ncop=min0(nadd,lentrackuse)
            incr0=max(ifirstframe-1,lastframe-ncop)
c            write (iout,9682) incr,incr0,nadd,ncop
c9682        format(' INCR=',I6,' INCR0=',I6,' NADD=',I6,' NCOP=',I6)
c           print *,'NCOP,NSEG_SCR=',NCOP,NSEG_SCR
            if (ncop .ge. nseg_scr .and. nseg_scr .gt. 0) then
              lenseg_scr=ncop/nseg_scr
              call scramble(ixran,nseg_scr)
              do is=1,nseg_scr-1
                ishift=(ixran(is)-is)*lenseg_scr
c               print *,'IS,IXRAN(IS),ISHIFT=',is,ixran(is),ishift
                do i=lenseg_scr*(is-1)+1,lenseg_scr*is
                  itrack(incr+i)=itrack(incr0+i+ishift)
                end do
              end do
c             Last (shorter) segment
              do i=lenseg_scr*(nseg_scr-1)+1,ncop
                itrack(incr+i)=itrack(incr0+i)
              end do
            else
              do i=1,ncop
                itrack(incr+i)=itrack(incr0+i)
              end do
            end if
            nadd=nadd-ncop
            incr=incr+ncop
          end do
        else if (iauctype .eq. 5) then
c         Pad with random states, p(on)=percent on
          nframeuse=lastframeinp-ifirstframe+1
          do i=lastframe+1,lastframeinp
            call randpx(1,ran)
            irand=0
            if (ran(1) .lt. percon) irand=1
            itrack(i)=irand
          end do
        end if
        lastframeuse=lastframeinp
      end if
c     print *,'LASTFRAMEUSE,IFR=',lastframeuse,ifr
      do ifr=ifirstframe,max0(lastframe,lastframeinp)
        if (itrack(ifr) .eq. 1) then
          ntodo=min0(nframeuse/2,lastframeuse-ifr+1)-1
          do ifrr=1,ntodo
            auc(ifrr)=auc(ifrr)+itrack(ifr+ifrr)
            nauc(ifrr)=nauc(ifrr)+1
c           if (auc(ifrr) .gt. nauc(ifrr)) write (6,7923)
c    -        ifr,ifrr,ntodo,auc(ifrr),nauc(ifrr),itrack(ifr+ifrr)
c7923       format(' ifr,ifrr=',2i6,' ntodo=',i6,' auc=',f12.1,
c    -    ' nauc=',i6,' itrack=',i9)
          end do
        end if
      end do
c     Normalize auc
      iauchalf=0
      aucsum=0.0
c      write (0006,9858) nframe2,nframeuse,lastframeuse
c9858  format(' NFRAME2=',i6,' NFRAMEUSE=',i6,' LASTFRAMEUSE=',I6)
      do ifr=1,nframe2
        if (nauc(ifr) .gt. 0) auc(ifr)=auc(ifr)/float(nauc(ifr))
c       if (auc(ifr) .lt. 0.0 .or. auc(ifr) .gt. 1.0) then
c         print *,'IFR=',ifr,' AUC=',auc(ifr),' NAUC=',nauc(ifr),
c    -     ' NFRAME2=',nframe2
c         stop
c       end if
        aucsum=aucsum+auc(ifr)
        if (iauchalf .eq. 0) then
          if (auc(ifr) .lt. 0.5) iauchalf=ifr
        end if
      end do
      rnfr2=framefac*(nframeuse/2)
      if (iframeunit .gt. 2) then
        if (iauchalf .eq. 0) write (iout,1002) itr,aucsum*framefac,'>',
     -    rnfr2,frunit(iframeunit)(1:lfrunit(iframeunit))
        if (iauchalf .gt. 0) write (iout,1002) itr,aucsum*framefac,'=',
     -    iauchalf*framefac,frunit(iframeunit)(1:lfrunit(iframeunit))
      else
        if (iauchalf .eq. 0) write (iout,1001) itr,aucsum*framefac,'>',
     -    rnfr2,frunit(iframeunit)(1:lfrunit(iframeunit))
        if (iauchalf .gt. 0) write (iout,1001) itr,aucsum*framefac,'=',
     -    iauchalf*framefac,frunit(iframeunit)(1:lfrunit(iframeunit))
      end if
      write (iout,1003) (auc(i),i=1,10)
      incr=max0(1,(nframe-1)/20+1)
      write (iout,1004) (auc(i*incr),i=1,nframeuse/(2*incr))
      if (itr_ix .le. MAXCOPY) then
        call trnsfr(auc_all(1,itr_ix),auc,nframe2)
        auc_all(nframe2+1,itr_ix)=aucsum
      else
        nauc_extra=nauc_extra+1
      end if
      return
1001  format(' BA#',i4,' AUC sum=',f7.1,' Half time of the track auc ',
     -  a,f9.1,1x,a)
1002  format(' BA#',i4,' AUC sum=',f11.4,' Half time of the track auc ',
     -  a,f9.4,1x,a)
1003  format(' AUC (1-10):',10f6.3)
1004  format(' AUC (incr):',10f6.3)
1005  format(' Paddig switched to 0s (nframe=',i6,' laston=',i6,
     -  ' loffmin=',i6,')')
      end
