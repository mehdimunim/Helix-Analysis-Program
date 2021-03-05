      subroutine clstrs_maxnn(ineig,nneig,n0,n,iclst,ifirst,
     -  ilast,nofcls0,nofcls,idrop,maxnode,maxgr,maxneig)
c#    MMC routine 462 lstmod: 07/15/02
c*****Find all clusters in a network and sort atoms in a cluster by groups
      dimension ineig(maxneig,maxnode),nneig(maxnode),iclst(maxnode),
     -  ifirst(maxgr),ilast(maxgr),idrop(maxnode)
      data indxidel /0/
c     Input parameters:
c     n0,n: Use atomnumbers (vertices) from n0 to n (inclusive)
c     maxnode,maxneig,maxgr: Array sizes - see dimension statement above
c     nneig(i) : number of neighbours (bonded) of node i
c     ineig(j,i) : j-th neighbour of atom i
c     Workspace array: idrop
c     Output parameters:
c     nofcls0: Number of disconnected clusters (groups) found previously
c     nofcls: Number of disconnected clusters (groups) found
c     iclst,ifirst,ilast: The elements of the ig-th group are atoms
c     (iclst(ia),ia=ifrst(ig),ilast(ig))
c
c     Description of the algorithm:
c     Pick the node with the largest # of neighbors (and lowest energy)
c     Make it a cluster; remove it ; repeat
c     Initialization
c     print *,'CLSTRS_nn start n0,n,maxgr=',n0,n,maxgr
      nleft=n-n0+1
      nofcls=nofcls0
      lastnode=n0-1
      do while (nleft .gt. 0)
c       Find node(s) with largest nn
        nnmax=0
        do i=n0,n
          if (nneig(i) .ge. nnmax) then
            nnmax=nneig(i)
            imin=i
          end if
        end do
        if (nnmax .lt. 0) write (6,1002) imin,nnmax
        nofcls=nofcls+1
        lastnode=lastnode+1
        ifirst(nofcls)=lastnode
        iclst(lastnode)=imin
        do in=1,nneig(imin)
         iclst(lastnode+in)=ineig(in,imin)
        end do
        lastnode=lastnode+nneig(imin)
        ilast(nofcls)=lastnode
        nleft=nleft-nneig(imin)-1
        if (ilast(nofcls)+nleft .ne. n)
     -    write (6,1001) nofcls,ilast(nofcls),nleft,n
c       write (6,1004) nofcls,ifirst(nofcls),ilast(nofcls),nleft,nnmax,
c    -    imin,nneig(imin)
        if (nleft .lt. 0) write (6,1000) nofcls0,nofcls,imin,nleft
c       Now eliminate imin and its neighbour from ineig
        call trnsfi(idrop,ineig(1,imin),nneig(imin))
        ndrop=nneig(imin)+1
        idrop(ndrop)=imin
        do i=1,ndrop
          id=idrop(i)
          if (nneig(id) .lt. 0) write (6,1003) nofcls,id,nneig(id)
          do in=1,nneig(id)
            idel=ineig(in,id)
            do jn=1,nneig(idel)
              if (ineig(jn,idel) .eq. id) then
                indxidel=jn
                go to 100
              end if
            end do
100         ineig(indxidel,idel)=ineig(nneig(idel),idel)
            nneig(idel)=nneig(idel)-1
          end do
          nneig(id)=-1
        end do
      end do
1000  format(' PROGRAM ERROR: nofcls0,nofcls,imin,nleft=',4i8,' < 0 !')
1001  format(' PROGRAM ERROR: nofcls,ilast(nofcls),nleft,n=',4i6,
     -  ' (il+nl ne n)')
1002  format(' PROGRAM ERROR: imin,nnmax=',2i5,' ( <0 !)')
1003  format(' PROGRAM ERROR: nofcls,id=',2i5,' nneig(id)=',i4,' (<0!)')
c1004  format(' nofcls,ifirst(nofcls),ilast(nofcls),nleft,nnmax=',5i5,/,
c     -  ' imin,nn(imin)=',2i5)
      return
      end
