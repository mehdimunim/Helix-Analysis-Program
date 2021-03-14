      subroutine checkwall(ioppbc,edge,cmin,cmax,c,cell,ncell1,ncell,
     -  ih,nslt,nnh,walld)
      dimension edge(3),cmin(3),cmax(3),c0(3),c(3,nslt),ih(nslt),
     -  cell(3,27)
      if (ioppbc .eq. 1 .or. ioppbc .eq. 2) then
        call extension(c,ih,nnh,1,nslt,cmin,cmax,c0,0,0,v)
        walld=1000.0
        do k=1,3
          wd=(edge(k)-(cmax(k)-cmin(k)))/2.0
          if (wd .lt. walld) walld=wd
        end do
      else
        walld=1.0
        do ia=1,nslt
          call genimdist(c(1,ia),cell,ncell1,ncell,icmin,rmin2)
          if (icmin .gt. 1) walld=-1.0
        end do
      end if
      return
      end
