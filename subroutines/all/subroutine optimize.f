      subroutine optimize(co,cn,cnew,ih,nslt,nnh,n,edge,edgen,edge_gen,
     -  ioppbc,npbc,ireduce,cell,ncell,cellalt,volnew,nwnew,
     -  closlim,closorg,closest,rot,ixyzhex,inn,jnn,ichk,iunsat,
     -  istickout,i2dopt,ftol)
      dimension co(3,n),cn(3,n),cnew(3,n),ih(n),edge(3),edgen(3),
     -  edge_gen(3,3),cell(3,27),cellalt(3,27),rot(3,3),ixyzhex(3)
c     Find the orientation maximizing image-image distances
      dimension cmin(3),cmax(3),edgenn(3)
      dimension simplex(4,3),angle(3),vertex(4)
c     print *,'OPTIMIZE start closlim,n,nslt=',closlim,n,nslt
      iunsat=0
      istickout=0
      call checkwall(ioppbc,edge,cmin,cmax,co,cell,1,ncell,
     -  ih,nslt,nnh,walldist)
      if (walldist .lt. 0.0) then
        volnew=1.e+33
        print *,'Optimization skipped because solute is out of the box'
        istickout=1
        return
      end if
c     print *,'edge=',edge,' ncell=',ncell
      closorg=distminimg(co,ih,nslt,nnh,edge,ioppbc,cell,ncell,ixyzhex,
     -  inno,jnno)
      if (nslt .lt. 3000 .and. ichk .eq. 1) then
        clo=distminck(co,ih,nslt,nnh,cell,cellalt,ioppbc,ncell,
     -    inno1,jnno1)
        if (abs(closorg-clo) .gt. 0.0001) print *,
     -    'PROGRAM ERROR: Original distance=',closorg,' Check=',clo
      end if
      if (i2dopt .eq. 0) then
c       3-D optimzation - initialize simplex
        call zeroit(simplex,12)
        do i=1,3
          simplex(i,i)=0.1
        end do
        do i=1,4
          do k=1,3
            angle(k)=simplex(i,k)
          end do
          vertex(i)=touch(co,ih,nslt,nnh,edge,ioppbc,cell,ncell,
     -      ixyzhex,angle,rot,1,0,cnew)
        end do
c       print *,'-Initial simplex distances=',vertex
c       call chkort(rot)
        call amoeba(simplex,vertex,4,3,3,ftol,iter,co,cnew,ih,nslt,nnh,
     -    edge,ioppbc,cell,ncell,ixyzhex,rot,1)
      else
c       2-D optimzation
        call opt2d(co,cnew,ih,n,nslt,nnh,edge,ioppbc,cell,ncell,
     -    ixyzhex,angle,rot,ftol,1,i2dopt,iter)
      end if
      print *,'--- Number of iterations=',iter
      call rotate_c(co,n,rot,cn,'OPTIMIZE',8)
      closest=distminimg(cn,ih,nslt,nnh,edge,ioppbc,cell,ncell,ixyzhex,
     -  inn,jnn)
      write (6,2006) closorg,closest
      if (closest .gt. closlim .and. ireduce .eq. 1) then
        if (npbc .gt. 1) then
          call cellreduce(0,cn,ih,nslt,nnh,edge,edgen,edge_gen,ioppbc,
     -      npbc,cell,ncell,cellalt,ixyzhex,closest,volnew,nwnew)
          call trnsfr(edgenn,edgen,3)
        else
          call trnsfr(edgenn,edge,3)
        end if
        call cellreduce(1,cn,ih,nslt,nnh,edgenn,edgen,edge_gen,ioppbc,
     -    npbc,cell,ncell,cellalt,ixyzhex,closlim,volnew,nwnew)
        write (6,2007) volnew,nwnew
c       Recreate the cell centers for edge
        call crorgn(edge,edge_gen,ioppbc,3,ncell,cell,cellalt,
     -    ixyzhex,rinscr,rcirc)
      else if (ireduce .eq. 1) then
        print *,'Smallest image distance is unsatisfactory:',closest
        iunsat=1
      end if
      return
2006  format(' Initial and optimized smallest image distances=',
     -  f5.2,f6.2,' A')
2007  format(' Smallest reduced volume=',f10.2,' A**3 accomodating ',i6,
     -  ' waters')
      end
