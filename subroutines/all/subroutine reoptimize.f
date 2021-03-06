      subroutine reoptimize(c,copt,c2,ih,nslt,nnh,n,edge,edgeopt,
     -  edge_gen,ioppbc,npbc,ireduce,cell,cellalt,ncell,volopt,volnew,
     -  closlim,closopt,closest,rot0,bestrot,inopt,jnopt,inn,jnn,
     -  ixyzhex,iunsat,istickout,impr,i2dopt,ftol)
      dimension c(3,n),copt(3,n),c2(3,n),ih(n),edge(3),edgeopt(3),
     -  edge_gen(3,3),edgen(3),cell(3,27),cellalt(3,27),rot0(3,3),
     -  bestrot(3,3),rotopt(3,3),ixyzhex(3)
c     print *,'Reoptimize n,ioppbc=',n,ioppbc
      impr=0
      volnew=volopt
      call optimize(c,c,c2,ih,nslt,nnh,n,edge,edgen,edge_gen,ioppbc,
     -  npbc,ireduce,cell,ncell,cellalt,volnew,nwnew,closlim,closorg,
     -  closest,rotopt,ixyzhex,inn,jnn,1,iunsat,istickout,i2dopt,ftol)
      if (istickout .gt. 0) return
      if (ireduce .eq. 1) then
c       Select smallest volume
        if (volnew .lt. volopt) then
          volopt=volnew
          closopt=closest
          inopt=inn
          jnopt=jnn
          call trnsfr(copt,c,3*n)
          call trnsfr(edgeopt,edgen,3)
          call matprod(rotopt,rot0,bestrot)
          impr=1
        end if
      else
c       Select largest minimum distance
        if (closest .gt. closopt) then
          closopt=closest
          inopt=inn
          jnopt=jnn
          call trnsfr(copt,c,3*n)
          call matprod(rotopt,rot0,bestrot)
          impr=1
        end if
      end if
      return
      end
