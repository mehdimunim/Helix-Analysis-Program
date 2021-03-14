      subroutine clstrs(ineig,nneig,nnloop,n0,n,iclst,ifirst,ilast,
     -  nofcls0,nofcls,iorder,iuseorder,iused,iclique,iout,inperr,iverb,
     -  maxat,maxgr,maxneig)
c#    MMC routine 462 lstmod: 07/15/02; clicque option added
c*****Find all clusters in a network and sort atoms in a cluster by groups
      dimension ineig(maxneig,maxat),nneig(maxat),iclst(maxat),
     -  ifirst(maxgr),ilast(maxgr),nnloop(maxat),iorder(maxat),
     -  iused(maxat)
c     Input parameters:
c     n0,n: Use atomnumbers (vertices) from n0 to n (inclusive)
c     maxat,maxneig,maxgr: Array sizes - see dimension statement above
c     nneig(i) : number of neighbours (bonded) of molecule i
c     ineig(j,i) : j-th neighbour of atom i
c     Workspace arrays:
c     iused(i) : 1 - atom i is not accounted for yet
c                0 - atom i is already accounted for
c     nnloop (i) : copy of nneig in loops 1 and 2 (temporary storage)
c                  the number of loop-closing bonds of atom i thereafter
c     Output parameters:
c     nofcls0: Number of disconnected clusters (groups) found previously
c     nofcls: Number of disconnected clusters (groups) found
c     iclst,ifirst,ilast: The elements of the ig-th group are atoms
c     (iclst(ia),ia=ifrst(ig),ilast(ig))
c
c     Description of the algorithm:
c     Starting with an atom, the algorithm successively includes its
c     neighbours and then the neighbours of the atoms already on the list.
c     By excluding atoms that 're-occurred' the algorithm essentially generates
c     a spanning tree.
c     Initialization
c     write (iout,*) 'CLSTRS n0,n,nofcls0,iverb=',n0,n,nofcls0,iverb
      jc=0
      do ii=n0,n
        if (iuseorder .eq. 1) then
          i=iorder(ii)
        else
          i=ii
        end if
        nnloop(i)=nneig(i)
        iused(i)=1
        if (iverb .gt. 0)
     -    write (iout,1010) i,nneig(i),(ineig(j,i),j=1,nneig(i))
        do in=1,nneig(i)
          j=ineig(in,i)
          if (j .lt. n0 .or. j .gt. n) then
            write (iout,1013) i,j,n0,n
            return
          end if
          nfound=0
          do jn=1,nneig(j)
            if (ineig(jn,j) .eq. i) nfound=nfound+1
          end do
          if (nfound .ne. 1) then
            write (iout,1011) j,i,i
            write (iout,1010) j,nneig(j),(ineig(k,j),k=1,nneig(j))
            inperr=inperr+1
          end if
        end do
      end do
      nofcls=nofcls0
      if (n .lt. n0) return
      nfound=n0-1
      ifirst(nofcls0+1)=n0
      do i=n0,n
        if (iused(i) .gt. 0) then
c         Start search
c         ncl is the number of elements in the cluster
c         ic is the index of the atom under consideration
c         kroot is the serial no of the lowest element in the cluster
c         that may still have neighbours not examined yet
          ic=i
          ncl=0
          kroot=1
c         NTOTSKIP=0
          do while (kroot .lt. ncl .or. ncl .eq. 0)
            if (iverb .gt. 1) write (iout,1009)
     -        i,kroot,ncl,ic,nnloop(ic),iused(ic)
            do while (ncl .eq. 0 .or. nnloop(ic) .gt. 0)
              if (iused(ic) .gt. 0) then
                if (iclique .eq. 0 .or. ncl .eq. nfound) then
c                 Include ic into the list
                  ncl=ncl+1
                  iused(ic)=0
                  iclst(nfound+ncl)=ic
                  if (iverb .gt. 1)
     -              write (iout,1008) 'Added1',ic,nnloop(ic),ncl,kroot
                else
c                 Check if ic is connected to all members found so far
                  nmatch=0
                  do icc=nfound+1,nfound+ncl
                    do iin=1,nneig(ic)
                      if (iclst(icc) .eq. ineig(iin,ic)) then
                        nmatch=nmatch+1
                        go to 101
                      end if
                    end do
c                   No match - skip
                    go to 102
101                 continue
                  end do
102               if (nmatch .eq. ncl) then
c                   New clique member found
                    ncl=ncl+1
                    iused(ic)=0
                    iclst(nfound+ncl)=ic
                    if (iverb .gt. 1)
     -                write (iout,1008) 'Added',ic,nnloop(ic),ncl,kroot
                  else
                    if (iverb .gt. 1) write (iout,1008)
     -                'Skipped',ic,nnloop(ic),ncl,kroot
                    iused(ic)=-1
c                   NTOTSKIP=NTOTSKIP+1
                  end if
                end if
c                write (iout,8734) (iclst(ii),ii=nfound+1,nfound+ncl)
c8734            format(' Clique so far:',(20i3))
              end if
              if (nnloop(ic) .gt. 0) then
c               Now search for the first unused neighbor of ic
                jc=ineig(nnloop(ic),ic)
                do while (iused(jc) .le. 0 .and. nnloop(ic) .gt. 0)
                  jc=ineig(nnloop(ic),ic)
                  nnloop(ic)=nnloop(ic)-1
                  if (iverb .gt. 1) write (iout,1006)
     -              ic,jc,nnloop(ic),iused(jc),ncl
                end do
                if (iused(jc) .gt. 0) ic=jc
              end if
c             if (NTOTSKIP .gt. 100) stop
            end do
c           Neighbour chain ended, back to the kroot-th atom in the list
            ic=iclst(nfound+kroot)
            if (nnloop(ic) .eq. 0) kroot=kroot+1
            if (iverb .gt. 1) write (iout,1007)
     -          nfound,kroot,ic,jc,nnloop(ic),ncl
          end do
c         Cluster of ncl elements found
c         stop
          nofcls=nofcls+1
          if (iclique .eq. 1) then
            do ii=n0,n
              nnloop(ii)=nneig(ii)
              if (iused(ii) .eq. -1) iused(ii)=1
            end do
          end if
c         memmax=0
c         memmin=10000000
c         do ia=nfound+1,nfound+ncl
c           if (memmax .lt. iclst(ia)) memmax=iclst(ia)
c           if (memmin .gt. iclst(ia)) memmin=iclst(ia)
c         end do
          nfound=nfound+ncl
          ilast(nofcls)=ifirst(nofcls)+ncl-1
          if (nofcls .lt. maxgr) ifirst(nofcls+1)=ilast(nofcls)+1
          if (iverb .gt. 0) write (iout,1005) nofcls,
     -      ifirst(nofcls),ilast(nofcls),ncl
        end if
      end do
      if (iverb .gt. 0) write (iout,*) 'Number of clusters=',nofcls
c     ilast(nofcls)=n
c     il=ifirst(nofcls)+ncl-1
      if (n .ne. ilast(nofcls)) then
        write (iout,1001) ilast(nofcls),n0,n
        inperr=inperr+1
      end if
      return
1001  format(' ***** PROGRAM ERROR in cluster search:',
     -  ' il=',i6,' n0,n=',2i6)
1005  format(' CLSTRS nofcls=',i8,' if,il=',2i6,' ncl=',i5)
1006  format(' CLSTRS ic,jc=',2i5,' nnloop(ic),iused(jc)=',2i3,
     -  ' ncl=',i5)
1007  format(' CLSTRS nfound,kroot=',2i6,' ic,jc,nnloop(ic)=',2i6,i3,
     -  ' ncl=',i5)
1008  format(' CLSTRS ',a,' ic=',i8,' nnloop(ic)=',i2,' ncl,kroot=',2i5)
1009  format(' CLSTRS i,kroot,ncl,ic,nnloop(ic),iused(ic)=',6i8)
1010  format(' CLSTRS ',i5,' nn=',i4,(' in=',15i6))
1011  format(' ***** PROGRAM ERROR: neighbour',i8,' of',i8,' does not ',
     -  'have',i8,' as neighbor exactly once:')
1013  format(' ERROR: index=',i5,' neighbor=',i5,' is outside the ',
     -  '[',i5,',',i5,'] interval')
      end
