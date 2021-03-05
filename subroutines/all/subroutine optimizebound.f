      subroutine optimizebound(co,cn,cnew,ih,nslt,nnh,n,rot,mintyp,
     -  ixyzhex,i2dopt,ftol)
      dimension co(3,n),cn(3,n),cnew(3,n),ih(n),rot(3,3)
c     Find the orientation minimizing the size of the bounding cube
      dimension edge(3),cell(3,27),ixyzhex(3)
      dimension simplex(4,3),angle(3),vertex(4)
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
     -      ixyzhex,angle,rot,mintyp,0,cnew)
        end do
c       print *,'-Initial simplex distances=',vertex
c       call chkort(rot)
        call amoeba(simplex,vertex,4,3,3,ftol,iter,co,cnew,ih,nslt,nnh,
     -    edge,ioppbc,cell,ncell,ixyzhex,rot,mintyp)
      else
c       2-D optimzation
        call opt2d(co,cnew,ih,n,nslt,nnh,edge,ioppbc,cell,ncell,
     -    ixyzhex,angle,rot,ftol,mintyp,i2dopt,iter)
      end if
      print *,'--- Number of iterations=',iter
      call rotate_c(co,n,rot,cn,'OPTBOUND',8)
      return
      end
