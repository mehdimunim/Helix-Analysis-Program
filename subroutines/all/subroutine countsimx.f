      subroutine countsimx(ifclst1,ilclst1,index1,nclust1,
     -  ifclst2,ilclst2,index2,nclust2,rmsdsim,iout,mx2d)
      dimension ifclst1(mx2d),ilclst1(mx2d),index1(mx2d),
     -  ifclst2(mx2d),ilclst2(mx2d),index2(mx2d)
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (4*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
     -  ihbpair(2,MAXBONDS),a1(MAXBONDS),a2(MAXBONDS),fill(IFILL2)
c     Count the  number of pairs in each cluster that are within rmsdsim
      write (iout,1000) nclust1,nclust2,rmsdsim
      do ic1=1,nclust1
        do ic2=1,nclust2
          nsimic=0
          rmsddmin=1000000.0
          rmsddmax=0.0
          do i1=ifclst1(ic1),ilclst1(ic1)
            do i2=ifclst2(ic2),ilclst2(ic2)
              rmsd=rmsd2d(index1(i1),index2(i2))
              if (rmsd .le. rmsdsim) nsimic=nsimic+1
              if (rmsd .lt. rmsddmin)rmsddmin=rmsd
              if (rmsd .gt. rmsddmax)rmsddmax=rmsd
            end do
          end do
          nmem1=ilclst1(ic1)-ifclst1(ic1)+1
          nmem2=ilclst2(ic2)-ifclst2(ic2)+1
          write (iout,1001) ic1,nmem1,ic2,nmem2
          if (nmem1*nmem2 .gt. 0 .and. min0(nmem1,nmem2) .gt. 1)
     -      write (iout,1002) rmsdsim,nsimic,
     -        100.0*float(nsimic)/float(nmem1*nmem2),rmsddmin,rmsddmax
        end do
        write (iout,*)
      end do
      return
1000  format(' Counting similar members between',i4,' and ',i4,
     -  ' clusters',/,' Similarity RMSD cutoff=',f5.1,' A')
1001  format(' Traj1, cluster',i4,' (nmem=',i4,') and ',
     -  ' Traj2, cluster',i4,' (nmem=',i4,')')
1002  format(3x,'# of Traj1-Traj2 pairs within',f6.1,' A RMSD=',i8,
     -  ' (',f6.2,' % of all pairs)',/,
     -  3x,'RMSD range of all Traj1-Traj2 pairs: [',f7.1,',',f7.1,']')
      end
