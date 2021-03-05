      subroutine clstrs_kmeans(n0,n,index2d,iclst,ifirst,ilast,nofcls,
     -  rmsd2d,icent,it1,it2,it3,it4,it5,indxclst,c,cent,cent_prev,
     -  ifail,maxnode,maxgr,mx2d)
c*****Find all clusters in a network and sort atoms in a cluster by groups
      dimension index2d(maxnode),iclst(maxnode),ifirst(maxgr),
     -  ilast(maxgr),it1(maxnode),it2(maxnode),it3(maxnode),
     -  it4(maxnode),it5(maxnode),indxclst(maxnode),c(3,n),
     -  cent(3,maxnode),cent_prev(3,maxnode)
      dimension icent(mx2d),rmsd2d(mx2d,mx2d)
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
c     k-means++
c     Cluster centers are the COM's of the cluster members
c     print *,'CLSTRS_K_COM start n0,n,nofcls=',n0,n,nofcls
c     Initialization
      iverb=0
      call indexit(index2d,1,MAX2D,0)
      do i=1,n
        do j=i+1,n
          rmsd2d(i,j)=dist2(c(1,j),c(1,i))
          rmsd2d(j,i)=rmsd2d(i,j)
        end do
      end do
      call init_kmedoids(nofcls,index2d,ifirst,ilast,icent,rmsd2d,n0,n,
     -  nnode,ierr,iverb,mx2d)
      ifail=0
      do i=1,nofcls
        call trnsfr(cent(1,i),c(1,icent(i)),3)
      end do
      call zeroiti(indxclst,n0-1,n)
      nchange=1
      iter=0
      maxiter=max0(5*nofcls,2*(n-n0))
      do while (nchange .gt. 0 .and. iter .lt. maxiter)
        nchange=0
        iter=iter+1
c       Find the nearest center to each node
        do i=n0,n
          rmin=10000.0
          jmin=0
          do j=1,nofcls
            d2=dist2(cent(1,j),c(1,i))
            if (d2 .lt. rmin) then
              rmin=d2
              jmin=j
            end if
          end do
          if (jmin .ne. indxclst(i)) then
            nchange=nchange+1
            indxclst(i)=jmin
          end if
        end do
c       write (6,7711) 'Bef sort indxclst:',(indxclst(i),i=n0,n)
c7711   format(' KMEANS ',a,(/i3,19i4))
        if (nchange .gt. 0) then
c         Find the new centers
          call trnsfr(cent_prev,cent,3*nofcls)
          call zeroit(cent,3*nofcls)
          call zeroiti(icent,0,nofcls)
          do i=n0,n
            ic=indxclst(i)
            icent(ic)=icent(ic)+1
            do k=1,3
              cent(k,ic)=cent(k,ic)+c(k,i)
            end do
          end do
          do ic=1,nofcls
            if (icent(ic) .gt. 0) then
              do k=1,3
                cent(k,ic)=cent(k,ic)/float(icent(ic))
              end do
            else
              call trnsfr(cent,cent_prev,3*nofcls)
            end if
          end do
        end if
      end do
      call indexit(it1,1,nnode,0)
c    do i=n0,n
c       w(i-n0+1)=indxclst(i)
c     end do
      call mrgsrti(6,it1,indxclst(n0),nnode,it2,it3,it4,it5,n)
c     write (6,7711) 'aft sort indxclst:',(it1(i),i=1,n-n0+1)
c     write (6,7711) 'aft sort it1:',(it1(i),i=1,n-n0+1)
      ifirst(1)=n0
      indxprev=indxclst(n0)
      nc=1
      do i=n0,n
        iclst(i)=it1(i-n0+1)
        if (indxclst(i) .gt. indxprev) then
          indxprev=indxclst(i)
          ilast(nc)=i-1
          nc=nc+1
          ifirst(nc)=i
        end if
      end do
      ilast(nc)=n
      if (nc .ne. nofcls) then
        write (6,1001) nc,nofcls
        ifail=1
      end if
      if (nchange .gt. 0) write (6,1002) iter
c     write (6,1000) iter,(icent(k),k=1,nofcls)
c      do ic=1,nofcls
c       if (ifirst(ic) .le. ilast(ic)) then
c        write (40,7721) ic,(iclst(ia),ia=ifirst(ic),ilast(ic))
c7721     format(' ic=',i3,/,(20i4))
c        else
c          write (40,*) 'Cluster ',ic,' is empty'
c      end do
      return
c1000  format(' Iteration ',i4,' Number of cluster members=',(10i5))
1001  format(' Clustering failure: nc=',i5,
     -  ' Number of clusters requested=',i5,/,' Clustering algorithm ',
     -  'needs extension',/' Try with different starting structure or ',
     -  'use the option where cluster centers are nodes')
1002  format(' NOTE: clustering did not converge in',i6,' iterations ')
      return
      end
