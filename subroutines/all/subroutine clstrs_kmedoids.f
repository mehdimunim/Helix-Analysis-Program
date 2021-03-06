      subroutine clstrs_kmedoids(n0,n,index2d,iclst,ifirst,ilast,nofcls,
     -  rmsd2d,icenttyp,icent,icent0,it1,it2,it3,it4,it5,indxclst,t1,
     -  icent_fin,maxnode,maxgr,mx2d,iout)
c*****Find all clusters in a network and sort atoms in a cluster by groups
      dimension index2d(maxnode),iclst(maxnode),ifirst(maxgr),
     -  ilast(maxgr),it1(maxnode),it2(maxnode),it3(maxnode),
     -  it4(maxnode),it5(maxnode),indxclst(maxnode),icent_fin(maxgr),
     -  icent(mx2d),icent0(mx2d),t1(maxnode),rmsd2d(mx2d,mx2d)
      parameter (MAX2D=5000)
      dimension indxmax(MAX2D),indxav(MAX2D),distmax(MAX2D),
     -  distav(MAX2D)
      character*1 ans
      data icminmax /0/,icminav /0/,jmin /0/
c     Input parameters:
c     n0,n: Use atomnumbers (vertices) from n0 to n (inclusive)
c     maxnode,maxgr: Array sizes - see dimension statement above
c     Workspace arrays: it1,it2,it3
c     Output parameters:
c     nofcls: Number of clusters requested
c     iclst,ifirst,ilast: The elements of the ig-th group are atoms
c     (iclst(ia),ia=ifrst(ig),ilast(ig))
c
c     Description of the algorithm:
c     k-means -> k-medoids
c     Cluster centers have to be one of the nodes.
c     The center node has the smallest max dev from the rest of the cluster
c     print *,'CLSTRS_K start n0,n,nofcls,mx2d=',n0,n,nofcls,mx2d
c     Initialization
c      write (6,9671) (index2d(i),i=n0,n)
c9671  format(' CLSTRS_KMEANS INDEX2D:',/,(20i4))
      iverb=1
      nitersave=min0(20,mx2d/nofcls)
      call zeroiti(icent,0,nitersave*nofcls)
      call init_kmedoids(nofcls,index2d,ifirst,ilast,icent,rmsd2d,n0,n,
     -  nnode,ierr,iverb,mx2d)
      call quiz(ans,icenttyp,'r',' ',0,'k-medoids center choice',23,0,
     -  5,6,139)
      if (icenttyp .eq. 1) write (iout,1004) 'smallest maximum distance'
      if (icenttyp .eq. 2) write (iout,1004) 'smallest average distance'
      call sortlist(0,icent,nofcls,it1,it2,'IK0',0,mx2d)
c     print *,'CENT=',(icent(i),i=1,nofcls)
      call zeroiti(indxclst,n0-1,n)
      do ic=1,nitersave
        icent0(ic)=(ic-1)*nofcls
      end do
      nchange=1
      iter=0
      maxiter=max0(5*nofcls,2*(n-n0))
      looping=0
      do while (nchange .gt. 0 .and.
     -          iter .lt. maxiter .and. looping .eq. 0)
        nchange=0
        ic0prev=icent0(mod(iter,nitersave)+1)
        iter=iter+1
        ic0=icent0(mod(iter,nitersave)+1)
c       Find the nearest center to each node
        do i=n0,n
           ii=index2d(i)
          jmin=0
          do j=1,nofcls
            if (i .eq. icent(ic0prev+j)) then
              jmin=j
            end if
          end do
          if (jmin .eq. 0) then
            rmin=10000.0
            do j=1,nofcls
              if (rmin .gt. rmsd2d(ii,index2d(icent(ic0prev+j)))) then
                rmin=rmsd2d(ii,index2d(icent(ic0prev+j)))
                jmin=j
              end if
            end do
          end if
          if (jmin .ne. indxclst(i)) then
            nchange=nchange+1
            indxclst(i)=jmin
          end if
        end do
c       write (6,7711) 'Bef sort indxclst:',(indxclst(i),i=n0,n)
c7711   format(' KMEANS ',a,(/i3,19i4))
        if (nchange .gt. 0) then
c         Find the new centers
          call indexit(it1,1,nnode,0)
c         do i=n0,n
c           w(i-n0+1)=indxclst(i)
c         end do
          call mrgsrti(6,it1,indxclst(n0),nnode,it2,it3,it4,it5,n)
c         write (6,7711) 'aft sort indxclst:',(it1(i),i=1,n-n0+1)
c         write (6,7711) 'aft sort it1:',(it1(i),i=1,n-n0+1)
          ifirst(1)=n0
          ixprev=indxclst(n0)
          nc=1
          do i=n0,n
            iclst(i)=it1(i-n0+1)
            if (indxclst(i) .gt. ixprev) then
              ixprev=indxclst(i)
              ilast(nc)=i-1
              nc=nc+1
              ifirst(nc)=i
            end if
          end do
          ilast(nc)=n
          nchange=0
          if (nc .ne. nofcls) then
            write (6,1001) nc,nofcls
          else
            do ic=1,nofcls
              rminmax=100000.0
              ravmax=100000.0
              ixa=0
              do ia=ifirst(ic),ilast(ic)
                iaa=index2d(iclst(ia))
                rmax=0.0
                rav=0.0
                do ja=ifirst(ic),ilast(ic)
                  if (rmsd2d(iaa,index2d(iclst(ja))) .gt. rmax)
     -              rmax=rmsd2d(iaa,index2d(iclst(ja)))
                    rav=rav+rmsd2d(iaa,index2d(iclst(ja)))
                end do
                if (rminmax .gt. rmax) then
                  rminmax=rmax
                  icminmax=ia
                end if
                rav=rav/float(ilast(ic)-ifirst(ic)+1)
                if (ravmax .gt. rav) then
                  ravmax=rav
                  icminav=ia
                end if
                ixa=ixa+1
                distmax(ixa)=rmax
                distav(ixa)=rav
                indxav(ixa)=ixa
                indxmax(ixa)=ixa
              end do
              nmem=ixa
              if (icenttyp .eq. 1) then
                icent(ic0+ic)=iclst(icminmax)
              else if (icenttyp .eq. 2) then
                icent(ic0+ic)=iclst(icminav)
              else if (icenttyp .eq. 3) then
                call mrgsrt(6,indxmax,distmax,nmem,it2,it3,it4,t1,nmem)
                call mrgsrt(6,indxav,distav,nmem,it2,it3,it4,t1,nmem)
c                write (77,9873) (indxmax(i),i=1,nmem)
c                write (77,9874) (distmax(i),i=1,nmem)
                do ia=1,nmem
                  it2(indxmax(ia))=ia
                  it3(indxav(ia))=ia
                end do
c                write (77,9875) (it2(i),i=1,nmem)
c                write (77,9876) (it3(i),i=1,nmem)
c9873            format(' INDXMAX:',/,(20i5))
c9874            format(' DISTMAX:',/,(20f5.1))
c9875            format(' IT2:',/,(20i5))
c9876            format(' IT3:',/,(20i5))
                icnext=0
                avminindx=2*MAX2D
                do ia=1,nmem
                  av=float(it2(ia)+it3(ia))/2.0
                  t1(ia)=av
                  if (av .lt. avminindx) then
                    avminindx=av
                    icnext=ia
                  end if
                end do
c                write (77,9877) (t1(i),i=1,nmem)
c9877            format(' RAVINDX:',/,(20f5.1))
c                write (77,*) 'ICNEXT=',icnext
                icent(ic0+ic)=iclst(ifirst(ic)-1+icnext)
              end if
            end do
            call mrgsortlist(icent(ic0+1),it1,it2,it3,it4,it5,nofcls)
            do ic=1,nnode
              if (icent(ic0prev+ic) .ne. icent(ic0+ic))nchange=nchange+1
            end do
            if (nchange .gt. 0 .and. iter .gt. nitersave) then
              ip=2
              do while (looping .eq. 0 .and. ip .lt. nitersave)
                idiff=0
                ic0prev=icent0(mod(iter-ip,nitersave)+1)
                ip=ip+1
                do ic=1,nofcls
                  if (icent(ic0prev+ic) .ne. icent(ic0+ic)) idiff=1
                end do
                looping=1-idiff
              end do
            end if
          end if
        else
          ic0=icent0(mod(iter-1,nitersave)+1)
        end if
        if (iverb .gt. 0) write (6,1000) iter,(icent(ic0+k),k=1,nofcls)
        do ic=1,nofcls
          icent_fin(ic)=icent(ic0+ic)
        end do
        write (iout,1005) (icent_fin(ic),ic=1,nofcls)
c       do ic=1,nofcls
c         write (iout,7721) ic,ifirst(ic),ilast(ic),
c    -      (iclst(ia),ia=ifirst(ic),ilast(ic))
c7721     format(' ic=',i3,' IFIRST,ILAST=',2i5,/,(20i4))
c       end do
      end do
      if (nchange .gt. 0 .and. looping .eq. 0) write (6,1002) iter
      if (nchange .gt. 0 .and. looping .eq. 1) write (6,1003) iter
      return
1000  format(' Iteration ',i4,' Centers chosen=',(10i5))
1001  format(' PROGRAM ERROR: nc=',i5,' Number of clusters requested=',
     -  i5)
1002  format(' NOTE: clustering did not converge in',i6,' iterations ',
     -  '- it may be looping')
1003  format(' NOTE: clustering detected looping on the centers found ',
     -  'after',i6,' iterations')
1004  format(' Criterion for cluster center choice: ',a)
1005  format(' Cluster centers:',(10i5))
      return
      end
