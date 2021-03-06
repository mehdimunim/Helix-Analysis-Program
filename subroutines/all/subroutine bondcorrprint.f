      subroutine bondcorrprint(nhbcorr,iout,line,index,iresno,indexbond,
     -  correxp,icorrtyp,icorrtrans,nframeav,nhbdist,ifhbclust,
     -  ilhbclust,ixhbclst,rclust,nclust,iclust,inamcol1,inamcol2,
     -  irescol1,irescol2,nhneigmin,hblimfac,angmin,rhph_sltbmax,
     -  percmin,percmax,minresdist,maxresdist,numres,icorrp,nframe,
     -  bondname,lbondname,icontract,ixresno,resnames,maxrec,mxbonds)
      dimension iresno(maxrec),index(maxrec),indexbond(mxbonds),
     -  nhbdist(mxbonds),ifhbclust(mxbonds),ilhbclust(mxbonds),
     -  ixhbclst(mxbonds),ixresno(maxrec)
      character*(*) bondname
      character*132 line(maxrec)
      character*8 resnames(maxrec)
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (6*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
     -  ihbtores(MAXBONDS),nusepair(MAXBONDS),nhb_atot(MAXBONDS),
     -  nhb_rtot(MAXBONDS),ibond(MAXBONDS),a1(MAXBONDS),fill(IFILL2)
      common /bondpairs/ ihbpair(2,MAXBONDS),ihb_pair_res(3,MAXBONDS)
      dimension i12(2),ihistogram(10)
      character*1 corb
      character*8 aname(2),rname(2)
c     print *,'BONDCORRP nclust,nhbcorr,iclust,bondname=',
c    -  nclust,nhbcorr,iclust,bondname(1:lbondname)
      if (nhbcorr .gt. MAX2D) then
        print *,'PROGRAM ERROR: nhbcorr=',nhbcorr,' MAX2D=',MAX2D
        return
      end if
      lan=inamcol2-inamcol1+1
      lrn=irescol2-irescol1+1
      if (iclust .eq. 0) write (iout,1026) bondname(1:lbondname)
      if (iclust .eq. 1) write (iout,1013) bondname(1:lbondname)
      if (icorrtyp .eq. 1) write (iout,1014)
      if (icorrtyp .eq. 2) write (iout,1015)
      if (icorrtyp .eq. 3) write (iout,1017) nframeav
      if (icorrtrans .eq. 1) write (iout,1018) correxp
      if (icorrtrans .eq. 2) write (iout,1019) correxp
      if (icorrtrans .eq. 3) write (iout,1020) correxp
      call write_traj_lim(iout,
     -  'Correlations calculated from trajectory',39,1,incr_tr,0)
      if (bondname(1:8) .eq. 'hydrogen') then
        write (iout,1010) hblimfac,angmin
      else if (bondname(1:11) .eq. 'hydrophobic') then
        write (iout,1011) bondname(1:11),rhph_sltbmax,' ',nhneigmin
      else if (bondname(1:11) .eq. 'salt-bridge') then
        write (iout,1011) bondname(1:11),rhph_sltbmax
      end if
      nprint=0
      if (percmin .gt. 0.0 .or. percmax .lt. 100.0) then
        write (iout,1006) bondname(1:lbondname),percmin,percmax
        nprint=1
      end if
      if (minresdist .gt. 1 .or. maxresdist .lt. numres-1) then
        write (iout,1007) bondname(1:lbondname),minresdist,maxresdist
        nprint=1
      end if
      if (nprint .eq. 0) write (iout,1005) bondname(1:lbondname)
      if (iclust .gt. 0) then
        if (rclust .gt. 0.0)  write (iout,1003) rclust
        nclp=nclust
      else
        nclp=1
        ifhbclust(1)=1
        ilhbclust(1)=nhbcorr
      end if
      maxhbtype=0
      minhbtype=nframe
      scpmax=0.0
      scpmin=1.0
      if (icorrtrans .ne. 2) then
        scpmaxp=1.0
      else
        scpmaxp=5.0
        if (correxp .gt. 1.0) scpmaxp=10.0
      end if
      call zeroiti(ihistogram,0,10)
      do ib=1,nhbcorr
        nhbd=nhbdist(indexbond(ib))
        if (nhbd .gt. maxhbtype) maxhbtype=nhbd
        if (nhbd .lt. minhbtype) minhbtype=nhbd
        do ibb=2,ib
          if (icontract .eq. 0) then
            scp=rmsd2d(indexbond(ib),indexbond(ibb))
          else
            scp=rmsd2d(ixhbclst(ib),ibb)
          end if
          if (scpmax .lt. scp) scpmax=scp
          if (scpmin .gt. scp) scpmin=scp
          ix=10*abs(scp-1.e-6)/scpmaxp+1
          if (ix .gt. 10) ix=1
          ihistogram(ix)=ihistogram(ix)+1
        end do
      end do
      write (6,1016) scpmin,scpmax
      write (iout,1016) scpmin,scpmax
      normfac=max0(1,(nhbcorr*(nhbcorr-1))/2)
      write (6,1008) scpmaxp,
     -   (float(ihistogram(i))/float(normfac),i=1,10)
      write (iout,1008) scpmaxp,
     -   (float(ihistogram(i))/float(normfac),i=1,10)
      ihb=0
      corb=':'
      if (icorrp .eq. 0) corb=' '
c      write (iout,9877) (indexbond(i),i=1,nhbcorr)
c9877  format(' INDEXBOND:',/,(20i3))
      do ic=1,nclp
        if (iclust .eq. 1) then
          write (iout,1012) ic,ilhbclust(ic)-ifhbclust(ic)+1
          call findbestcorrep(iout,ifhbclust(ic),ilhbclust(ic),
     -      indexbond,mxbonds)
        end if
        if (icontract .eq. 0) write (iout,1022) bondname(1:lbondname)
        if (icontract .eq. 1) write (iout,1023)
        if (icorrp .gt. 0) write (iout,1024)
        do ib=ifhbclust(ic),ilhbclust(ic)
          ihb=ihb+1
          if (icorrp .eq. 1) write (iout,*)
          ib1=indexbond(ib)
          if (icontract .eq. 0) then
            do i=1,2
              i12(i)=ihbpair(i,ib1)
              aname(i)(1:lan)=line(index(i12(i)))(inamcol1:inamcol2)
              rname(i)(1:lrn)=line(index(i12(i)))(irescol1:irescol2)
            end do
            write (iout,1001) ib,ib1,(i12(i),
     -        aname(i)(1:lan),iresno(i12(i)),rname(i)(1:lrn),i=1,2),
     -        float(100*nhbdist(ib1))/float(nframe),corb
          else
            ixr1=ihb_pair_res(1,ib1)
            ixr2=ihb_pair_res(2,ib1)
            ir1=ixresno(ixr1)
            ir2=ixresno(ixr2)
            call lastchar(resnames(ixr1),lc1,8)
            call lastchar(resnames(ixr2),lc2,8)
            write (iout,1025) ixhbclst(ib),ib1,resnames(ixr1)(1:lc1),
     -        ir1,resnames(ixr2)(1:lc2),ir2
          end if
          if (icorrp .eq. 1) then
            if (icontract .eq. 0) then
              if (icorrtrans .eq. 2) then
                write (iout,1021)
     -            (rmsd2d(ib1,indexbond(ib2)),ib2=1,nhbcorr)
              else
                write (iout,1009)
     -            (rmsd2d(ib1,indexbond(ib2)),ib2=1,nhbcorr)
              end if
            else
              if (icorrtrans .eq. 2) then
                write (iout,1021)
     -            (rmsd2d(ixhbclst(ib),ib2),ib2=1,nhbcorr)
              else
                write (iout,1009)
     -            (rmsd2d(ixhbclst(ib),ib2),ib2=1,nhbcorr)
              end if
              if (icorrtrans .eq. 1) then
                do ib2=1,nhbcorr
                 if (rmsd2d(ixhbclst(ib),ib2) .lt. 0.0 .or.
     -               rmsd2d(ixhbclst(ib),ib2) .gt. 1.0) write (iout,*)
     -             'ib2=',ib2,' rmsd2d=',rmsd2d(ixhbclst(ib),ib2)
                end do
              end if
            end if
          end if
        end do
        if (nclp .gt. 1) write (iout,1002) ic
      end do
      return
1001  format(' FB#',i5,' (',i4,'):',i6,1x,a,i6,1x,a,' -',
     -  i6,1x,a,i6,1x,a,' occ=',f6.2,'% ',a)
1002  format(1x,23('-'),' End of cluster ',i4,1x,30('-'))
1003  format(' Clustered with Rclust=',f5.3)

1005  format(' All ',a,' bonds are included')
1006  format(' Filtered list: ',a,' bonds are included only when ',/,
     -  5x,'their occurrence is in the range [',f5.1,' - ',f6.1,'] %')
1007  format(' Filtered list: ',a,' bonds are included only when ',/,
     -  5x,'their interresidue number difference is in the range [',i5,
     -  ' - ',i5,'] residues')
1008  format(' Histogram of distance measures (10% bins in the 0 - ',
     -  f4.1,' range):',/,10f6.3,/)
1009  format(15f5.2)
1010  format(' Hydrogen-bonds thresholds: hblimfac=',f5.2,
     -  ' and H-bond angle minimum=',f5.1)
1011  format(' Distance threshold for ',a,'bond=',f5.2,a,/,' Minimum ',
     -  'number of bonded hydrogens for a hydrophobic carbon=',i1)
1012  format(' Cluster ',i4,' number of members:',i5)
1013  format(' Clustering results of ',a,' bonds based on correlation')
1014  format(' Correlation measure: CORR=(HB1.HB2)+(1-HB1).(1-HB2) ',
     -  '(include the off states)')
1015  format(' Correlation measure: CORR=(HB1.HB2) (correlate only the',
     -  '  off states')
1016  format(' Range of the distance measure values:',
     -  ' [',f5.2,' - ',f8.2,']')
1017  format(' Correlations are calculated using a running average ',
     -  'over',i4,' windows')
1018  format(' Distance measure: the correlation complement power ',
     -  '(1-CORR)^(',f5.3,')')
1019  format(' Distance measure: the correlation inverse power ',
     -  '(1/CORR)^(',f5.3,')')
1020  format(' Distance measure: the correlation power CORR^(',f5.3,')')
1021  format(15f5.1)
1022  format(1x,a,' bond# (original #)')
1023  format(' res-res# (original #)')
1024  format(' Correlation measures in each row refer to the ',
     -  'unclustered order')
1025  format(' RR#',i5,' (',i5,') ( ',a,i4,' - ',a,i4,')')
1026  format(' Correlation based distance matrix of ',a,' bonds')
      end
