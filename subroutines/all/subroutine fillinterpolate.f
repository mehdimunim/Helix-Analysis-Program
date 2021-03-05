      subroutine fillinterpolate(pmap,nx1,nx,ny1,ny,nz1,nz,igincr,
     -  nxmax,nymax,nzmax,pinf,ninterpol)
      dimension pmap(nxmax,nymax,nzmax)
c     print *,'FI nx1,nx,ny1,ny,nz1,nz=',nx1,nx,ny1,ny,nz1,nz
      ninterpol=0
      ngr=0
      do iz=nz1,nz,igincr
        do iy=ny1,ny,igincr
          do ix=nx1,nx,igincr
            ngr=ngr+1
            if (pmap(ix,iy,iz) .eq. pinf) then
              ninterpol=ninterpol+1
              pint=0.0
              ixlim=ix
              do while (ixlim .le. nx .and. pmap(ixlim,iy,iz) .eq. pinf)
                ixlim=ixlim+1
              end do
              if (ix .eq. nx1) then
                if (ixlim .gt. nx) then
c                 Empty line - stop for now
                  print *,'Empty line iy=',iy,' iz=',iz
                  stop
                end if
                pint=pint+pmap(ixlim,iy,iz)
              else if (ixlim .gt. nx) then
                pint=pint+pmap(ix-igincr,iy,iz)
              else
                pint=pint+(igincr*pmap(ixlim,iy,iz)+
     -            (ixlim-ix)*pmap(ix-igincr,iy,iz))/(ixlim-ix+igincr)
              end if
              iylim=iy
              do while (iylim .le. ny .and. pmap(ix,iylim,iz) .eq. pinf)
                iylim=iylim+1
              end do
              if (iy .eq. ny1) then
                if (iylim .gt. ny) then
c                 Empty line - stop for now
                  print *,'Empty line ix=',ix,' iz=',iz
                  stop
                end if
                pint=pint+pmap(ix,iylim,iz)
              else if (iylim .gt. ny) then
                pint=pint+pmap(ix,iy-igincr,iz)
              else
                pint=pint+(igincr*pmap(ix,iylim,iz)+
     -            (iylim-iy)*pmap(ix,iy-igincr,iz))/(iylim-iy+igincr)
              end if
              izlim=iz
              do while (izlim .le. nz .and. pmap(ix,iy,izlim) .eq. pinf)
                izlim=izlim+1
              end do
              if (iz .eq. nz1) then
                if (izlim .gt. nz) then
c                 Empty line - stop for now
                  print *,'Empty line ix=',ix,' iy=',iy
                  stop
                end if
                pint=pint+pmap(ix,iy,izlim)
              else if (izlim .gt. nz) then
                pint=pint+pmap(ix,iy,iz-igincr)
              else
                pint=pint+(igincr*pmap(ix,iy,izlim)+
     -            (izlim-iz)*pmap(ix,iy,iz-igincr))/(izlim-iz+igincr)
              end if
              pmap(ix,iy,iz)=pint/3.0
            end if
          end do
        end do
      end do
      return
      end
