      subroutine rmsdsubcluster(numclst,nfrst,n,index2d,iclst,ifirst,
     -  ilast,rmsdclust,isubclustertyp,surfacefract,denexphalf,
     -  i1,i2,i3,i4,i5,value,t1,iout,nofcls,label,llabel)
      dimension index2d(n),iclst(n),ifirst(n),ilast(n),value(n),t1(n),
     -  i1(n),i2(n),i3(n),i4(n),i5(n)
      character*(*) label
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (4*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
     -  ihbpair(2,MAXBONDS),a1(MAXBONDS),a2(MAXBONDS),fill(IFILL2)
      dimension ix1(MAX2D)
      data indxidel /0/
c     Break up a cluster by removing low-neighbor number nodes
      write (6,1002) nfrst,n,rmsdclust
      write (iout,1002) nfrst,n,label(1:llabel),rmsdclust
c     First recreate the nn list for the cluster selected in a way
c     to stay within the cluster
      nofcls=0
      call zeroiti(nng,0,n)
      do ii=nfrst,n
        i=index2d(ii)
        do jj=ii+1,n
          j=index2d(jj)
          if (rmsd2d(i,j) .lt. rmsdclust) then
            nng(ii)=nng(ii)+1
            nng(jj)=nng(jj)+1
            ing(nng(ii),ii)=jj
            ing(nng(jj),jj)=ii
          end if
        end do
      end do
      call clstrs(ing,nng,ix1,nfrst,n,iclst,ifirst,ilast,0,nofcls,i2,0,
     -  i2,0,6,inperr,0,n,n,MAX2D)
      if (nofcls .gt. 1) then
        write (6,1004) nofcls
        return
      end if
      call trnsfi(i3,nng,n)
      do while (.true.)
        istart=0
        if (isubclustertyp .eq. 1) then
c         Subclustering by # of neighbors
          call trnsfi(i4,nng,n)
c         Find NN range (singletons excluded)
          nnmax=-1
          nnmin=n+1
          do i=nfrst,n
            if (nng(i) .gt. 0) then
              if (nng(i) .gt. nnmax) nnmax=nng(i)
              if (nng(i) .lt. nnmin) nnmin=nng(i)
            end if
          end do
          if (nnmin .eq. nnmax) then
            write (6,1001) nnmin
            go to 200
          end if
c         Now eliminate all nodes with nnmin neighbors
          write (6,1007) nnmin
          write (iout,1007) nnmin
          ndel=0
          do i=nfrst,n
            if (nng(i) .eq. nnmin) then
              ndel=ndel+1
              i5(ndel)=i
            end if
          end do
          write (6,1003) numclst,n-nfrst+1,nnmin,nnmax
          write (iout,1003) numclst,n-nfrst+1,nnmin,nnmax
        else
c         Subclustering by 'density'
          tmin=100000.0
          tmax=0.0
          do ii=nfrst,n
            i=index2d(ii)
            t1(ii)=0.0
            do jj=nfrst,n
              j=index2d(jj)
              if (i .ne. j) t1(ii)=t1(ii)+1.0/rmsd2d(i,j)**denexphalf
            end do
            if (tmin .gt. t1(ii)) tmin=t1(ii)
            if (tmax .lt. t1(ii)) tmax=t1(ii)
          end do
          ndel=0
          do i=nfrst,n
            if (t1(i) .lt. surfacefract*tmax) then
              ndel=ndel+1
              i5(ndel)=i
            end if
          end do
          write (6,1009) numclst,n-nfrst+1,tmin,tmax
          write (iout,1009) numclst,n-nfrst+1,tmin,tmax
        end if
c       i5 has the list of nodes to delete
        write (6,1010) (index2d(i5(i)),i=1,ndel)
        write (iout,1010) (index2d(i5(i)),i=1,ndel)
        do i=1,ndel
          id=i5(i)
          do in=1,nng(id)
            idel=ing(in,id)
            do jn=1,nng(idel)
              if (ing(jn,idel) .eq. id) then
                indxidel=jn
                go to 100
              end if
            end do
100         ing(indxidel,idel)=ing(nng(idel),idel)
            nng(idel)=nng(idel)-1
          end do
          nng(id)=0
        end do
        call clstrs(ing,nng,i1,nfrst,n,iclst,ifirst,ilast,0,nofcls,
     -    i2,0,i2,0,6,inperr,0,n,n,MAX2D)
c       Count the number of non-singleton clusters
        nsingl=0
        do ic=1,nofcls
          if (ifirst(ic) .eq. ilast(ic)) then
            nsingl=nsingl+1
          end if
        end do
        nofcls=nofcls-nsingl
        write (6,1005) nofcls
        nc=0
        if (nofcls .gt. 1) then
          do ic=1,nofcls+nsingl
            if (ifirst(ic) .lt. ilast(ic)) then
              nc=nc+1
              i5(nc)=ilast(ic)-ifirst(ic)+1
            end if
          end do
        end if
        if (nc .gt. 0) write (6,1006) (i5(i),i=1,nc)
        call askyn('Do you want to shave more nodes',31,
     -    1,-1,moreshave,0,0)
        if (moreshave .eq. 0) go to 200
      end do
200   if (nofcls .gt. 1) then
c       Add back to the real clusters the deleted nodes
        call zeroiti(i5,0,n)
        do ic=1,nofcls
          if (ifirst(ic) .lt. ilast(ic)) then
            do ia=ifirst(ic),ilast(ic)
             i5(ia)=ic
            end do
          end if
        end do
c       Now i5 is the cluster number of a remaining node or zero
        do ia=nfrst,n
          if (i5(ia) .eq. 0) then
c           Find a remaining neighbor, assign ia to its cluster
            do ja=1,i3(ia)
              if (i5(ing(ja,i3(ia))) .ne. 0) then
                i5(ia)=i5(ing(ja,i3(ia)))
                go to 300
              end if
            end do
300         continue
          end if
        end do
c       Now, recluster the nodes based on i5
        call indexit(i1,1,n,0)
        nmem=n-nfrst+1
        do ia=1,nmem
          value(ia)=i5(ia+nfrst-1)
        end do
        call mrgsrt(6,i1,value,nmem,i2,i3,i5,t1,nmem)
        ic=1
        ifirst(ic)=nfrst
        icprev=value(1)
        do ia=nfrst,n
          ia0=ia-nfrst+1
          iclst(ia)=i1(ia0)
          if (value(ia0) .ne. icprev) then
            ilast(ic)=ia-1
            ic=ic+1
            ifirst(ic)=ia
          end if
        end do
        ilast(ic)=n
        if (ic .ne. nofcls) write (6,1000) nofcls,ic
c       Finally, rearrange the elements of index2d in the order of iclst
        do ia=nfrst,n
          i1(ia)=index2d(iclst(ia))
        end do
        call trnsfi(index2d(nfrst),i1(nfrst),n-nfrst+1)
      end if
      return
1000  format(' PROGRAM ERROR in rmsdsubcluster nofcls=',i5,' ic=',i5)
1001   format(' All nodes have the same number of neighbors (',i5,
     -  ') - can not subcluster')
1002  format(' Subclustering structures [',i6,' - ',i6,
     -  '] with ',a,' threshold=',f5.2)
1003  format(' Cluster #',i4,' (',i4,' members): minimum and maximum ',
     -  '# of neighbors=',2i5)
1004  format(' PROGRAM ERROR: invalid subclustering - nofcls=',i4)
1005  format(' Number of subclusters found=',i4)
1006  format(' Number of members in the subclusters found=',5i5,/(10i5))
1007  format(' Temporarily deleting members with only',i3,' neighbours')
1009  format(' Cluster #',i4,' (',i4,' members): local density ',
     -  'min/max=',2f10.5)
1010  format(' Temporarily deleting:',10i5)
      end
