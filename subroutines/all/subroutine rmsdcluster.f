      subroutine rmsdcluster(rmsdclust,nfrst,n,index2d,iwt,iclst,ifirst,
     -  ilast,it1,it2,it3,it4,it5,it6,c,cent,cent_prev,icent_fin,
     -  maxct,iclstyp,idenclstyp,nnmin,nofcls,nosortindex2d,iuseindex2d,
     -  ifail,ntry,label,llabel,iout)
      dimension index2d(n),iwt(n),iclst(n),ifirst(n),ilast(n),it1(n),
     -  it2(n),it3(n),it4(n),it5(n),it6(n),c(3,n),cent(3,maxct),
     -  cent_prev(3,maxct),icent_fin(maxct)
      character*(*) label
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL21=MAXPHI*MAXPHI*MAXPHI-
     -  (6*MAXBONDS+2*MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
     -  ihbtores(MAXBONDS),nusepair(MAXBONDS),nhb_atot(MAXBONDS),
     -  nhb_rtot(MAXBONDS),t1(MAXBONDS),w(MAXBONDS),nng2r(MAX2D),
     -  fill(IFILL21)
      dimension icent(MAX2D),icent0(MAX2D),index2dd(MAX2D),
     -  iorder(MAXBONDS),t2(MAXBONDS)
c     Set up neighbor list based on the threshold
c     print *,'RMSDCLUSTER iout,n,iclstyp,MX2D=',iout,n,iclstyp,MAX2D
c     print *,'RMSDCLUSTER iclstyp,idenclstyp,nnmin=',iclstyp,idenclstyp,nnmin
c     write (40,*)'RMSDCLUSTER  nfrst,n,iclstyp=',nfrst,n,iclstyp
c     write (40,*)'RMSDCLUSTER  RMSDCLUST=',rmsdclust
c      write (40,7934) (index2d(i),i=1,n)
c7934  format(' RMSDCLUSTER INDEX2D=',15i4)
      ifail=0
      sdmin=100000.0
      sdmax=0.0
      nz_sd=0
      nofcls_t=0
      if (nfrst .gt. n) then
        print *,'PROGRAM ERROR: nfrst > n - nfrst=',nfrst,' n=',n
        stop
      end if
      if (n .gt. MAXBONDS) then
        write (6,1001)
        stop
      end if
      if (nfrst .eq. n) then
c        Cluster has one element
         nofcls=1
         iclst(n)=n
         ifirst(1)=n
         ilast(1)=n
         return
      end if
      call zeroiti(nng,0,n)
      if (iclstyp .eq. 9) call zeroiti(nng2r,0,n)
      if (iuseindex2d .eq. 1) then
        call trnsfi(index2dd,index2d,n)
      else
        call indexit(index2dd,1,n,0)
      end if
      rmsdclust4=rmsdclust*4.0
c     write (77,*) 'RMSDCLUSTER nfrst,n=',nfrst,n
c     write (79,*) 'nfrst,n=',nfrst,n
c     do i=1,n
c       write (79,1004) (rmsd2d(i,j),j=1,n)
c     end do
c1004 format(5e13.6)
      do ii=nfrst,n
        i=index2dd(ii)
        do jj=ii+1,n
          j=index2dd(jj)
          if (rmsd2d(i,j) .lt. sdmin) sdmin=rmsd2d(i,j)
          if (rmsd2d(i,j) .gt. sdmax) sdmax=rmsd2d(i,j)
          if (rmsd2d(i,j) .eq. 0.0) nz_sd=nz_sd+1
c         if (rmsd2d(i,j) .eq. 0.0) then
c           write (6,*) 'ii,jj=',ii,jj,'i,j=',i,j
c         end if
          if (iclstyp .ne. 4 .and. iclstyp .ne. 5) then
            if (rmsd2d(i,j) .lt. rmsdclust) then
              nng(ii)=nng(ii)+1
              nng(jj)=nng(jj)+1
              ing(nng(ii),ii)=jj
              ing(nng(jj),jj)=ii
            end if
          end if
          if (iclstyp .eq. 9) then
            if (rmsd2d(i,j) .lt. rmsdclust4) then
              ing(MAX2D-nng2r(ii),ii)=jj
              ing(MAX2D-nng2r(jj),jj)=ii
              nng2r(ii)=nng2r(ii)+1
              nng2r(jj)=nng2r(jj)+1
            end if
          end if
        end do
      end do
      write (6,1002) sdmin,sdmax
      write (iout,1002) sdmin,sdmax
      if (iclstyp .eq. 8 .and. nosortindex2d .eq. 0) then
c       Clique-clustering - sort lists
c       Sort iorder by iwt
        nmem=n-nfrst+1
        call indexit(iorder,1,n,0)
        do ii=nfrst,n
          t2(ii-nfrst+1)=-float(iwt(ii))
        end do
        call mrgsrt(6,iorder(nfrst),t2,nmem,it2,it3,it4,t1,nmem)
c       Sort neighbors by iwt
        do i=nfrst,n
          if (nng(i) .gt. 1) then
            do ii=1,nng(i)
              t2(ii)=-float(iwt(ing(ii,i)))
            end do
            call mrgsrt(6,ing(1,i),t2,nng(i),it2,it3,it4,t1,nng(i))
          end if
        end do
      end if
c      do ii=nfrst,n
c        i=index2dd(ii)
c        write (40,8967) ii,i,(rmsd2d(i,index2dd(j)),j=1,n)
c8967    format(' ii,i=',2i5,' RMSD2D:',15f5.2,/,(15f5.2))
c        write (40,8968) ii,i,(ing(jj,ii),jj=1,nng(ii))
c8968    format(' ii,i=',2i5,' ING:',15i4,/,(14i4))
c      end do
c      do i=1,n
c        write (6,8701) i,nng(i),(ing(j,i),j=1,nng(i))
c8701    format(i4,' nn=',i4,' in=',10i5)
c      end do
      if (sdmin .eq. 0.0) then
        if (ntry .lt. 2)
     -    write (6,1000) label(1:llabel),'minimum','some',' ',nz_sd
      else if (sdmax .eq. 0.0) then
        write (6,1000) label(1:llabel),'maximum','all'
        print *,'Exiting clustering'
        ifail=1
        return
      end if
      if (iclstyp .le. 3) then
c       Single-link clustering
        call clstrs(ing,nng,it1,nfrst,n,iclst,ifirst,ilast,0,nofcls,
     -    it2,0,it2,0,iout,inperr,0,n,n,MAX2D)
      else if (iclstyp .eq. 4) then
c       K-medoids clustering - no COM
        call clstrs_kmedoids(nfrst,n,index2dd,iclst,ifirst,ilast,nofcls,
     -    rmsd2d,icenttyp,icent,icent0,it1,it2,it3,it4,it5,it6,t1,
     -    icent_fin,n,n,MAX2D,iout)
        if (maxct .gt. 1) then
c         Extract coordinates of center nodes
          do ic=1,nofcls
            call trnsfr(cent(1,ic),c(1,icent_fin(ic)),3)
          end do
        end if
      else if (iclstyp .eq. 5) then
c       K-means clustering -  COM-based
        call clstrs_kmeans(nfrst,n,index2dd,iclst,ifirst,ilast,nofcls,
     -    rmsd2d,icent,it1,it2,it3,it4,it5,it6,c,cent,cent_prev,ifail,
     -    n,n,MAX2D)
        if (ifail .gt. 0) return
      else if (iclstyp .eq. 6 .or. iclstyp .eq. 7) then
        call clstrs_maxnn(ing,nng,nfrst,n,iclst,ifirst,ilast,0,
     -    nofcls,it2,n,n,MAX2D)
      else if (iclstyp .eq. 8) then
c       Clique-based clustering
        iverb=1
        iout_c=iout
c       iout_c=6
        call clstrs(ing,nng,it1,nfrst,n,iclst,ifirst,ilast,0,nofcls,
     -    iorder,1,it2,1,iout_c,inperr,iverb,n,n,MAX2D)
      else if (iclstyp .eq. 9) then
        call clstrs_density(nfrst,n,iclst,ifirst,ilast,nofcls,nofcls_t,
     -    nng,ing,nng2r,idenclstyp,nnmin,it1,it2,it3,it4,it5,it6,MAX2D,
     -    iout)
      end if
c     Sort iclslt within each cluster
c      write (06,9671) 'U',nofcls,(iclst(i),i=1,n)
c9671  format(1x,a,' NCLUST=',i4,' ICLST:',/,(20i4))
      do ic=1,nofcls
        nmem=ilast(ic)-ifirst(ic)
        if (nmem .ge. 1) then
          do ia=1,nmem
            t2(ia)=iclst(ia+ifirst(ic)-1)
          end do
          call indexit(it1,1,nmem,0)
          call mrgsrt(6,it1,t2,nmem,it2,it3,it4,t1,nmem)
          do ia=1,nmem
            iclst(ia+ifirst(ic)-1)=t2(ia)
          end do
        end if
      end do
c     Rearrange the elements of index2d in the order of iclst
      do ii=nfrst,n
        if (iclst(ii) .lt. 1 .or. iclst(ii) .gt. n)
     -    print *,'ii=',ii,' iclst=',iclst(ii)
        it1(ii)=index2d(iclst(ii))
      end do
      call trnsfi(index2d(nfrst),it1(nfrst),n-nfrst+1)
      if (iclstyp .eq. 9) nofcls=nofcls_t
      return
1000  format(1x,a,1x,a,' is zero - ',a,' items are identical',a,/,
     -  ' # of zero RMSDs=',i6)
1001  format(' ERROR: number of items exceeds limit',/,
     -  ' Redimension the arrays index2dd and iorder to M A X R E C')
1002  format(' Range of the distance measures: [',f10.5,',',f10.5,']')
      end
