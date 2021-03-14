      subroutine bondcorrsum(nbondcorr,nresbondcorr,scpmin,scpmax,
     -  ibframe,itempsum,ifirstframe,correxp,icorrtyp,
     -  icorrtrans,iaggregate,icorrbothstart,nframeav,imatprint,iout,
     -  mx2d)
      dimension ibframe(mx2d),itempsum(mx2d),ifirstframe(mx2d)
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
      character*1 ans
c     print *,'BONDCORRSUM nbondcorr,nframe,nresbondcorr=',
c    -  nbondcorr,nframe,nresbondcorr
      if (nbondcorr .le. 1) then
        write (6,1004) nbondcorr
        return
      end if
      ncorr=nbondcorr
      if (iaggregate .gt. 0) ncorr=nresbondcorr
      if (ncorr .gt. MAX2D) then
        write (6,1005) MAX2D,ncorr
        return
      end if
      call trajlimtest(nframe,MAXFRAMES)
      call quiz(ans,icorrtyp,'b',' ',0,'bond correlation type',21,
     -  0,5,6,85)
      if (icorrtyp .eq. 3) then
8010    call getint('Number of frames to average over',32,5,1,nframe,
     -    nframeav,84)
        if (nframeav .ge. nframe) then
          write (6,1008) nframe
          go to 8010
        end if
      end if
      call quiz(ans,icorrtrans,'-',' ',0,
     -  'correlation - measure transformation',36,0,5,6,000)
      write (6,1002)
      call getint('Exponent of correlations',24,1,0,9,icorrexp,3)
      if (nframe .ge. 1000)
     -  print *,'Calculating correlation matrix - wait'
      correxp=icorrexp
      if (icorrexp .lt. 1) correxp=1.0/float(-icorrexp)
      write (iout,1001)
      call zeroiti(nng,0,MAX2D)
      call zeroiti(ifirstframe,0,MAX2D)
      if (icorrtyp .eq. 3) call zeroiti(itempsum,0,MAX2D)
      call zeroiti(ing,0,MAX2D*MAX2D)
      nframe10=nframe/10+1
      scpmax=0.0
      scpmin=1.0
c     if (iaggregate .gt. 0) then
c       ihbtoresmax=0
c       do ib=1,nbondcorr
c         if (ihbtores(ib) .gt. ihbtoresmax) ihbtoresmax=ihbtores(ib)
c       end do
c     end if
      do ifr=1,nframe
        call readbitc(ires(1,ifr),ibframe,nbondcorr,30,MAXITEMS)
        if (icorrtyp .eq. 3) then
          do ib=1,ncorr
            itempsum(ib)=itempsum(ib)+ibframe(ib)
          end do
        end if
        do ib=1,ncorr
          if (ibframe(ib) .eq. 1 .and. ifirstframe(ib) .eq. 0)
     -      ifirstframe(ib)=ifr
          nng(ib)=nng(ib)+ibframe(ib)
        end do
c       Accumulate scalar product
        do ib1=2,ncorr
          do ib2=1,ib1-1
            if (icorrbothstart .eq. 0) then
               iuseib12=ifirstframe(ib1)+ifirstframe(ib2)
            else
               iuseib12=ifirstframe(ib1)*ifirstframe(ib2)
            end if
            if (iuseib12 .gt. 0) then
              if (icorrtyp .lt. 3) then
                icorr=ibframe(ib1)*ibframe(ib2)
                if (icorrtyp .eq. 1) then
c                 Both on and off states are correlated
                  icorr=icorr+(1-ibframe(ib1))*(1-ibframe(ib2))
                end if
                ing(ib1,ib2)=ing(ib1,ib2)+icorr
                ing(ib2,ib1)=ing(ib2,ib1)+1
              else if (mod(ifr,nframeav) .eq. 0) then
                icorr=itempsum(ib1)*itempsum(ib2)
                ing(ib1,ib2)=ing(ib1,ib2)+icorr
                ing(ib2,ib1)=ing(ib2,ib1)+1
              end if
            end if
          end do
        end do
        if (icorrtyp .eq. 3) then
          if (mod(ifr,nframeav) .eq. 0) call zeroiti(itempsum,0,ncorr)
        end if
        if (nframe .ge. 1000 .and. mod(ifr,nframe10) .eq. 0)
     -    write (6,1000) ifr/nframe10
      end do
c     Prepare the distance measure matrix
      scp=0.0
      do ib1=2,ncorr
        do ib2=1,ib1-1
          if (icorrtyp .eq. 1 .or. icorrtyp .eq. 2) then
            scp=float(ing(ib1,ib2))/float(ing(ib2,ib1))
          else if (icorrtyp .eq. 3) then
c           Average states
            scp=float(ing(ib1,ib2))/float(ing(ib2,ib1)*nframeav**2)
          end if
          if (scp .lt. 0.0 .or. scp .gt. 1.0)
     -      write (iout,1003) ib1,ib2,ing(ib1,ib2),ing(ib2,ib1),scp
          if (icorrtrans .eq. 1) then
            scp=(1.0-scp)**correxp
          else if (icorrtrans .eq. 2) then
            scp=(amin1(1.0/scp,99.0))**correxp
          else
            scp=scp**correxp
          end if
c          if (iaggregate .eq. 1)
c     -      write (iout,8781) ib1,ib2,ing(ib1,ib2),ing(ib2,ib1),scp
c8781      format(' IB1,2=',2i6,' ing(ib1,ib2)=',i6,
c     -      ' ing(ib2,ib1)=',i6,' scp=',f10.5)
          rmsd2d(ib1,ib2)=scp
          rmsd2d(ib2,ib1)=scp
          if (scp .lt. scpmin) scpmin=scp
          if (scp .gt. scpmax) scpmax=scp
        end do
      end do
      if (imatprint .eq. 1) then
        if (iaggregate .eq. 0) write (iout,1007) ' '
        if (iaggregate .eq. 1) write (iout,1007) 'residue-aggregated'
        do ib1=1,ncorr
          write (iout,1006) ib1,(rmsd2d(ib1,ib2),ib2=1,ncorr)
        end do
      end if
c      write (iout,8791) (nng(i),i=1,ncorr)
c      write (iout,8792) (ifirstframe(i),i=1,ncorr)
c8791  format(' NNG:',10i6)
c8792  format(' IFF:',10i6)
      return
1000  format(i2,'0 % done')
1001  format(/,' === Calculating correlation between bond histories ',
     -  'for clustering the bonds')
1002  format(' Negative number -n for exponent will be changed to ',
     -  '1/n')
1003  format(' PROGRAM ERROR: ib1,2=',2i6,' ing(ib1,ib2)=',i6,
     -  ' ing(ib2,ib1)=',i6,' scp=',f10.5)
1004  format(' Correlation calculations require at least 2 bonds',/,
     -  ' but there is only ',i2,' left')
1005  format(' ERROR: Bond correlation calculation is limited to ',i5,
     -  ' bonds',/,8x,'- filter out more bonds or recompile Simulaid ',
     -  'with MAX2D > ',i5)
1006    format(' IB1=',i4,' RMSD2D:',/,(15f5.2))
1007  format(' The ',a,' bond correlation matrix')
1008  format(' # of frames to average over should be < # of frames (',
     -  i6,')')
      end
