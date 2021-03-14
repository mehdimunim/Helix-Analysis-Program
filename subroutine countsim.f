      subroutine countsim(ifclst,ilclst,index,nclust,rmsdsimdefr,
     -  rmsdsim,nsim,iout,mx2d)
      dimension ifclst(mx2d),ilclst(mx2d),index(mx2d),nsim(mx2d)
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (4*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
     -  ihbpair(2,MAXBONDS),a1(MAXBONDS),a2(MAXBONDS),fill(IFILL2)
c     Count the  number of pairs in each cluster that are within rmsdsim
      rmsdsimdef=rmsdsimdefr
      if (rmsdsimdefr .eq. 0.0) rmsdsimdef=1.0
      call getreal('MAXimum RMSD for intra-cluster similarity',41,
     -  rmsdsimdef,rmsdsim,1,86)
      write (iout,1000) nclust,rmsdsim
      do ic=1,nclust
        nsimic=0
        rmsdmin=99999.9
        rmsdmax=0.0
        do ia=ifclst(ic),ilclst(ic)
          do ja=ia+1,ilclst(ic)
            rmsd=rmsd2d(index(ia),index(ja))
            if (rmsd .le. rmsdsim) nsimic=nsimic+1
            if (rmsd .lt. rmsdmin) rmsdmin=rmsd
            if (rmsd .gt. rmsdmax) rmsdmax=rmsd
          end do
        end do
        nsim(ic)=nsimic
        nmem=ilclst(ic)-ifclst(ic)+1
        write (iout,1001) ic,nmem
        if (nmem .gt. 1)
     -    write (iout,1002) rmsdsim,nsimic,
     -      100.0*float(nsimic)/float((nmem*(nmem-1))/2),rmsdmin,rmsdmax
      end do
      return
1000  format(/,' Counting similar members of',i4,' clusters with ',
     -  'similarity RMSD cutoff=',f5.1,' A')
1001  format(' Cluster',i5,' nmem=',i4)
1002  format(8x,'# of member pairs within',f6.1,' A RMSD=',i8,
     -  ' (',f6.2,' % of all pairs)',/,
     -  8x,'RMSD range of all pairs: [',f7.1,',',f7.1,']')
      end
