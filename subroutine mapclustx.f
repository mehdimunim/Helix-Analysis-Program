      subroutine mapclustx(ixtr2,ifclst1,ilclst1,irepmx1,nclust1,
     -  ifclst2,ilclst2,index2,nclust2,nndist,nndistsum,
     -  rmsdmapmax,nframe2,trajnam1,ltrajnam1,trajnam2,
     -  ltrajnam2,iout,mx2d)
      dimension ifclst1(mx2d),ilclst1(mx2d),irepmx1(mx2d),ifclst2(mx2d),
     -  ilclst2(mx2d),index2(mx2d),nndist(mx2d),nndistsum(mx2d)
      character*(*) trajnam1,trajnam2
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (4*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
     -  ihbpair(2,MAXBONDS),a1(MAXBONDS),a2(MAXBONDS),fill(IFILL2)
      dimension rmsddmin(MAX2D),rmsddmax(MAX2D)
      data iclmin /0/
c     Map traj 2 ont clusters of traj1
c     print *,'MAPCLUSTX nclust1,nclust2,iout=',nclust1,nclust2,iout
      write (iout,1004) rmsdmapmax,trajnam1(1:ltrajnam1),
     -  trajnam2(1:ltrajnam2)
      call zeroiti(nndistsum,0,nclust1)
      do ic2=1,nclust2
        do ic1=1,nclust1
          nndist(ic1)=0
          rmsddmax(ic1)=0.0
          rmsddmin(ic1)=1000000.0
        end do
        do ia=ifclst2(ic2),ilclst2(ic2)
          rmsddmn=100000.0
          if (ixtr2 .eq. 2) then
            do ic1=1,nclust1
              if (rmsd2d(irepmx1(ic1),index2(ia)) .lt. rmsddmn) then
                rmsddmn=rmsd2d(irepmx1(ic1),index2(ia))
                iclmin=ic1
              end if
            end do
          else
            do ic1=1,nclust1
              if (rmsd2d(index2(ia),irepmx1(ic1)) .lt. rmsddmn) then
                rmsddmn=rmsd2d(index2(ia),irepmx1(ic1))
                iclmin=ic1
              end if
            end do
          end if
          if (rmsddmn .le. rmsdmapmax)then
            nndist(iclmin)=nndist(iclmin)+1
            if (rmsddmn .lt. rmsddmin(iclmin)) rmsddmin(iclmin)=rmsddmn
            if (rmsddmn .gt. rmsddmax(iclmin)) rmsddmax(iclmin)=rmsddmn
          end if
        end do
        write (iout,1000) ic2,ilclst2(ic2)-ifclst2(ic2)+1,ixtr2,3-ixtr2
        do ic1=1,nclust1
          if (nndist(ic1) .gt. 0) then
            write (iout,1001) nndist(ic1),
     -        float(100*nndist(ic1))/float(ilclst2(ic2)-ifclst2(ic2)+1),
     -        ic1,ilclst1(ic1)-ifclst1(ic1)+1,
     -        ' ',rmsddmin(ic1),rmsddmax(ic1)
            nndistsum(ic1)=nndistsum(ic1)+nndist(ic1)
          else
            write (iout,1001) nndist(ic1),0.0,ic1,
     -        ilclst1(ic1)-ifclst1(ic1)+1
          end if
          write (iout,1005) ic2,irepmx1(ic2),ic1,irepmx1(ic1),
     -      rmsd2d(irepmx1(ic1),irepmx1(ic2))
        end do
      end do
      write (iout,*)
      nmapped=0
      do ic1=1,nclust1
        write (iout,1002) 3-ixtr2,ic1,ixtr2,nndistsum(ic1),
     -    100.0*float(nndistsum(ic1))/float(nframe2)
        nmapped=nmapped+nndistsum(ic1)
      end do
      write (iout,1003) ixtr2,3-ixtr2,nmapped,
     -  100.0*float(nmapped)/float(nframe2)
      return
1000  format(' Mapping cluster',i3,' (# of members=',i4,') of Traj',i2,
     -  ' to clusters of Traj',i2)
1001  format(i6,' members (',f5.1,'%) mapped on Clst',i4,
     -  ' (',i5,' members)',a,'RMSDs: [',f5.1,',',f5.1,']')
1002  format(' Traj',i2,' Clst',i4,': total # of ',
     -  'structures in Traj',i2,' mapped on it:',i5,' (',f5.1,' %)')
1003  format(' Total number of structures in Traj',i2,' mapped on Traj',
     -  i2,':',i4,' (',f5.1,' %)')
1004  format(' Mapping clusters of two trajectories:',/,
     -  ' counting the number of members of cluster A that are within',
     -  f6.2,' A RMSD',/,' of the representative element of cluster B',
     -  /,' Trajectory 1:',a,/,' Trajectory 2:',a)
1005  format(' RMSD between representative element of cluster',i3,' (',
     -  i6,') and',/,' representative element of cluster',i3,' (',i6,
     -  '):',f5.1,' A')
      end
