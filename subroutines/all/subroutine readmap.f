      subroutine readmap(iu,xstart,ystart,zstart,gx,gy,gz,ngx,ngy,ngz,
     -  nconf,work,lwork)
      dimension work(lwork)
      parameter (MAXPHI=400)
      common /nnwork/ phimap(MAXPHI,MAXPHI,MAXPHI)
      character*60 toplbl !ascii header
      integer*4 ivary ! 0 => x index varys most rapidly
      integer*4 nbyte ! =4, # of bytes in data
      integer*4 inddat ! =0, floating point data
      real*4 xang,yang,zang ! =90,90,90 unit cell angles
      integer*4 intx,inty,intz ! =igrid-1, # of intervals/grid side
      real*4 extentx, extenty, extentz ! maximum extent of grid in x,y,z
      real*4 xstart,xend ! beginning, end of grid sides
      real*4 ystart,yend ! in fractional
      real*4 zstart,zend ! units of extent*
      data imin /0/,jmin /0/,kmin /0/
      rewind iu
      read (iu) toplbl
c     print *,'READING MAP ',toplbl
      read (iu) ivary, nbyte, inddat, extentx, extenty, extentz,
     -  xang, yang, zang, xstart, xend, ystart, yend,
     -  zstart, zend, intx, inty, intz
      if (ivary .eq. 0) write (6,*) 'X index varies faster (Fortran)'
      if (ivary .gt. 0) write (6,*) 'Z index varies faster (C)'
c     print *,'ivary, nbyte, intdat=',ivary, nbyte, intdat
      gx=extentx*(xend-xstart)/(intx)
      gy=extenty*(yend-ystart)/(inty)
      gz=extentz*(zend-zstart)/(intz)
      phimin=1000000.0
      ngx=intx+1
      ngy=inty+1
      ngz=intz+1
      write (6,1002) 'X',extentx*xstart,extentx*xend,gx,ngx
      write (6,1002) 'Y',extenty*ystart,extenty*yend,gy,ngy
      write (6,1002) 'Z',extentz*zstart,extentz*zend,gz,ngz
      write (6,1003) extentx*(xstart+xend)/2.0,
     -  extenty*(ystart+yend)/2.0,extentz*(zstart+zend)/2.0
      nskip=0
      if (ngx .gt. MAXPHI) then
        ngxx=ngx
        do while (ngxx .gt. MAXPHI)
          ngxx=ngxx/2
          nskip=nskip+1
        end do
        if (nconf .le. 1) write (6,1000) ngx,MAXPHI,2**nskip
        gx=gx*2**nskip
        gy=gy*2**nskip
        gz=gz*2**nskip
      end if
      nuse=2**nskip
      if (ivary .eq. 0) then
        do k=1,ngz
          do j=1,ngy
            read (iu) (work(i),i=1,ngx)
            do i=1,ngx
              if (work(i) .lt. phimin) then
                phimin=work(i)
                imin=i
                jmin=j
                kmin=k
              end if
            end do
            if (mod(k,nuse) .eq. 0 .and. mod(j,nuse) .eq. 0) then
              do i=1,ngx
                if (mod(i,nuse) .eq. 0)
     -            phimap((i-1)/nuse+1,(j-1)/nuse+1,(k-1)/nuse+1)=work(i)
              end do
            end if
          end do
        end do
      else
        do k=1,ngx
          do j=1,ngy
            read (iu) (work(i),i=1,ngz)
            do i=1,ngz
              if (work(i) .lt. phimin) then
                phimin=work(i)
                imin=i
                jmin=j
                kmin=k
              end if
            end do
            if (mod(k,nuse) .eq. 0 .and. mod(j,nuse) .eq. 0) then
              do i=1,ngz
                if (mod(i,nuse) .eq. 0)
     -            phimap((i-1)/nuse+1,(j-1)/nuse+1,(k-1)/nuse+1)=work(i)
              end do
            end if
          end do
        end do
      end if
      xstart=xstart*extentx
      ystart=ystart*extenty
      zstart=zstart*extentz
      write (6,1001) phimin,xstart+(imin-1)*gx/2**nskip,
     -  ystart+(jmin-1)*gy/2**nskip,zstart+(kmin-1)*gz/2**nskip
      return
1000  format(' The number of gridpoints (',i4,') exceeds ',i4,
     -  ' - only every',i3,'-th will be used',/,
     -  ' (update/replace the nnwork common block to avoid this)')
1001  format(' The lowest potential=',f10.4,' at (',f8.2,',',f8.2,',',
     -  f8.2,')')
1002  format(1x,a,'-direction start=',f10.5,' end=',f10.5,
     -  ' gridsize=',f6.3,' (',i3,' grids)')
1003  format(' The grid is centered at (',f10.5,',',f10.5,',',f10.5,')')
      end
