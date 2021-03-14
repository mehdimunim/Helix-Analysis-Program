      subroutine delphigrid(iu,iu1,iu2,c,n,nslt,xstart,ystart,zstart,
     -  gx,gy,gz,ngx,ngy,ngz,igincr,rnear,igridfile,iexcl,iquery,
     -  interpol)
      dimension c(3,n)
c     Print the potential at the grid points
      parameter (MAXPHI=400)
      common /nnwork/ phimap(MAXPHI,MAXPHI,MAXPHI)
      character*1 xyz
      common /axislab/ xyz(3)
      dimension xyzmin(3),xyzmax(3)
      if (iquery .eq. 1) then
        x=0.0
        do while (x .ne. 999.0)
          call getreal('X coordinate (999 to quit)',26,999999.0,x,0,0)
          if (x .ne. 999.0) then
            call getreal('Y coordinate',12,999999.0,y,0,0)
            call getreal('Z coordinate',12,999999.0,z,0,0)
            call interpolate(x,y,z,gx,gy,gz,xstart,ystart,zstart,phi)
          write (6,1000) x,y,z,phi
          end if
        end do
      end if
      if (igridfile .eq. 1) then
c       Establish grid range to print
101     call getreal('Grid X coordinate minimum',25,xstart,xgmin,0,0)
        call getreal('Grid X coordinate maximum',25,xstart+ngx*gx,xgmax,
     -    0,0)
        ixgmin=(xgmin-xstart)/gx-1
        if (ixgmin .lt. 1) ixgmin=1
        ixgmax=(xgmax-xstart)/gx+1
        if (ixgmax .gt. ngx) ixgmax=ngx
        if (ixgmin .gt. ixgmax) then
          print *,'Invalid X grid range:',ixgmin,ixgmax
          go to 101
        end if
102     call getreal('Grid Y coordinate minimum',25,ystart,ygmin,0,0)
        call getreal('Grid Y coordinate maximum',25,ystart+ngy*gy,ygmax,
     -    0,0)
        iygmin=(ygmin-ystart)/gy-1
        if (iygmin .lt. 1) iygmin=1
        iygmax=(ygmax-ystart)/gx+1
        if (iygmax .gt. ngy) iygmax=ngy
        if (iygmin .gt. iygmax) then
          print *,'Invalid Y grid range:',iygmin,iygmax
          go to 102
        end if
103     call getreal('Grid Z coordinate minimum',25,zstart,zgmin,0,0)
        call getreal('Grid Z coordinate maximum',25,zstart+ngz*gz,zgmax,
     -    0,0)
        izgmin=(zgmin-zstart)/gz-1
        if (izgmin .lt. 1) izgmin=1
        izgmax=(zgmax-zstart)/gx+1
        if (izgmax .gt. ngz) izgmax=ngz
        if (izgmin .gt. izgmax) then
          print *,'Invalid Z grid range:',izgmin,izgmax
          go to 103
        end if
        rnear2=0.0
        if (iexcl .eq. 1) then
          rnear2=rnear**2
          ncheck=0
          ndrop=0
          do k=1,3
            xyzmin(k)=10000.0
            xyzmax(k)=-xyzmin(k)
          end do
          ii=0
          do ia=1,nslt
            do k=1,3
              if (c(k,ia) .lt. xyzmin(k)) xyzmin(k)=c(k,ia)
              if (c(k,ia) .gt. xyzmax(k)) xyzmax(k)=c(k,ia)
            end do
c           Mark grids to be dropped
            if (c(1,ia)+rnear .gt. xgmin .and.
     -          c(1,ia)-rnear .lt. xgmax .and.
     -          c(2,ia)+rnear .gt. ygmin .and.
     -          c(2,ia)-rnear .lt. ygmax .and.
     -          c(3,ia)+rnear .gt. zgmin .and.
     -          c(3,ia)-rnear .lt. zgmax) then
              ncheck=ncheck+1
              ixmin=(c(1,ia)-rnear-xstart)/gx-1
              if (ixmin .lt. 1) ixmin=1
              ixmax=(c(1,ia)+rnear-xstart)/gx+2
              if (ixmax .gt. ngx) ixmax=ngx
              iymin=(c(2,ia)-rnear-ystart)/gy-1
              if (iymin .lt. 1) iymin=1
              iymax=(c(2,ia)+rnear-ystart)/gy+2
              if (iymax .gt. ngy) iymax=ngy
              izmin=(c(3,ia)-rnear-zstart)/gz-1
              if (izmin .lt. 1) izmin=1
              izmax=(c(3,ia)+rnear-zstart)/gz+2
              if (izmax .gt. ngz) izmax=ngz
              do iz=izmin,izmax
                do iy=iymin,iymax
                  do ix=ixmin,ixmax
                    ii=ii+1
                    x=xstart+(ix-1)*gx
                    y=ystart+(iy-1)*gy
                    z=zstart+(iz-1)*gz
                    dd=(x-c(1,ia))**2+(y-c(2,ia))**2+(z-c(3,ia))**2
                    if (dd. lt. rnear2 .and.
     -                  phimap(ix,iy,iz) .ne. 999999.9) then
                      phimap(ix,iy,iz)=999999.9
                      ndrop=ndrop+1
                    end if
c                  write (77,7878) ia,ix,iy,iz,x,y,z,(c(k,ia),k=1,3),
c     -              dd,rnear2,ndrop
c7878              format(' ia=',i5,2x,3i5,' xyz=',3f10.4,' c=',3f10.4,
c     -              ' dd,rnear2=',2f10.2,' ndrop=',i9)
                  end do
                end do
              end do
            end if
          end do
          write (6,2008) (xyz(k),xyzmin(k),xyzmax(k),k=1,3)
        end if
        ngw=0
        ngd=0
        write (6,1002) ncheck,nslt,ii,ndrop
        do iz=izgmin,izgmax,igincr
          do iy=iygmin,iygmax,igincr
            do ix=ixgmin,ixgmax,igincr
              ngw=ngw+1
              x=xstart+(ix-1)*gx
              y=ystart+(iy-1)*gy
              z=zstart+(iz-1)*gz
              pm=phimap(ix,iy,iz)
              if (pm .eq. 999999.9) then
                pm=0.0
                ngd=ngd+1
              end if
              if (abs(pm) .lt. 99.0) then
                write (iu,1004) iz,x,y,z,pm
              else
                if (pm .lt. -999.0) pm=-999.0
                if (pm .gt. 999.0) pm=999.0
                write (iu,1010) iz,x,y,z,pm
              end if
              if (pm .ne. 0.0 .and. iu1 .gt. 0) then
                if (abs(pm) .lt. 99.0) then
                  write (iu1,1004) iz,x,y,z,pm
                else
                  if (pm .lt. -999.0) pm=-999.0
                  if (pm .gt. 999.0) pm=999.0
                  write (iu1,1010) iz,x,y,z,pm
                end if
              end if
            end do
          end do
        end do
        if (interpol .eq. 1) then
          call fillinterpolate(phimap,ixgmin,ixgmax,iygmin,iygmax,
     -     izgmin,izgmax,igincr,MAXPHI,MAXPHI,MAXPHI,999999.9,ninterpol)
          do iz=izgmin,izgmax,igincr
            do iy=iygmin,iygmax,igincr
              do ix=ixgmin,ixgmax,igincr
                x=xstart+(ix-1)*gx
                y=ystart+(iy-1)*gy
                z=zstart+(iz-1)*gz
                pm=phimap(ix,iy,iz)
                if (abs(pm) .lt. 99.0) then
                  write (iu2,1004) iz,x,y,z,pm
                else
                  if (pm .lt. -999.0) pm=-999.0
                  if (pm .gt. 999.0) pm=999.0
                  write (iu2,1010) iz,x,y,z,pm
                end if
              end do
            end do
          end do
          write (6,1005) 'interpolated',ninterpol
          if (ninterpol .ne. ngd)
     -      print *,'ERROR in number of dropped/interpolated grids'
        end if
      end if
      write (6,1003) ngw
      if (iexcl .gt. 0)
     -   write (6,1005) 'dropped (or written with zero potentials)',ngd
      return
1000  format(' phi(',f10.5,',',f10.5,',',f10.5,')=',f12.5)
1002  format(' Checked ',i5,' atoms out of ',i6,'; Checked',i9,
     -  ' grids',/,' Number of grids found too close to atoms=',i8)
1003  format(' Number of grid points written=',i7)
1004  format('ATOM        He   GRD A ',i3,4x,3f8.3,'  1.00',f6.2)
1005  format(' Number of ',a,' grid points=',i6)
1010  format('ATOM        He   GRD A ',i3,4x,3f8.3,'  1.00',f6.1)
2008  format(' Solute extensions: ',/,3(1x,a,': [',f7.1,' - ',f7.1,']'))
      return
      end
