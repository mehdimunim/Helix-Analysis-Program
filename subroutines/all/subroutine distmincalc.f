      subroutine distmincalc(ioppbc,cell,ncell,ixyzhex,edge,dc1,dc2,
     -  dc3,icell,sum)
      dimension cell(3,ncell),ixyzhex(3),edge(3)
      common /distmindata/ ef1,ef2,ef3,ef4,edgep,edgep2,edgey,edgex,
     -  edgex2,w,h,et1,et2,et3,et4,edge2(3),edghalf,edgsq3p2,edgp2sq3,
     -  edgsq2p3,edgpsq3,edg2psq3
      common /numbers/ sq3,sq3inv,sq3p2,sq2p3
      dimension dc123(3),d(3),summ(7)
      if (ioppbc .eq. 3) then
c       Face-centered cubic
        sum=dc1**2+dc2**2+dc3**2
        icell=1
        d1=abs(dc1)
        if (abs(dc2) .gt. d1) d1=abs(dc2)
        if (abs(dc3) .gt. d1) d1=abs(dc3)
        d1=-d1*ef2+ef4
        d2=abs(dc1)+abs(dc2)
        if (abs(dc2)+abs(dc3) .gt. d2)
     -    d2=abs(dc2)+abs(dc3)
        if (abs(dc1)+abs(dc3) .gt. d2)
     -    d2=abs(dc1)+abs(dc3)
        d2=-d2*ef1+ef3
        if (d1 .lt. d2) then
          sum1=sum+d1
        else
          sum1=sum+d2
        end if
        if (sum .gt. sum1) then
          sum=sum1
          icell=2
        end if
      else if (ioppbc .eq. 4 .or. ioppbc .eq. 5) then
c       Hexagonal prizm, prism along ixyzhex(1) ax, vertex on ixyzhex(2) ax
        dc123(1)=dc1
        dc123(2)=dc2
        dc123(3)=dc3
        dcc1=abs(dc123(ixyzhex(1)))
        dcc2=abs(dc123(ixyzhex(2)))
        dcc3=abs(dc123(ixyzhex(3)))
        icell=2
        if (dcc3 .gt. edgex2) then
c         Cell 6 or right side of cell 5
          if (h+sq3inv*(dcc3-edgex2) .le. dcc2) then
c           Cell 5
            dcc3=dcc3-edgex2
            dcc2=dcc2-w
          else
c           Cell 6
            dcc3=dcc3-edgex
          end if
        else
c         Cell 1 or left side of cell 5
          if (h+sq3inv*(edgex2-dcc3) .le. dcc2) then
c           Cell 5
            dcc3=dcc3-edgex2
            dcc2=dcc2-w
          else
            icell=1
          end if
        end if
        if (dcc1 .gt. edgep2) then
          dcc1=dcc1-edgep
          icell=2
        end if
        sum=dcc1**2+dcc2**2+dcc3**2
      else if (ioppbc .eq. 1 .or. ioppbc .eq. 2) then
c       Rectangular
        dcc1=abs(dc1)
        dcc2=abs(dc2)
        dcc3=abs(dc3)
        icell=1
        if (dcc1 .gt. edge2(1)) then
          dcc1=dcc1-edge(1)
          icell=2
        end if
        if (dcc2 .gt. edge2(2)) then
          dcc2=dcc2-edge(2)
          icell=2
        end if
        if (dcc3 .gt. edge2(3)) then
          dcc3=dcc3-edge(3)
          icell=2
        end if
        sum=dcc1**2+dcc2**2+dcc3**2
      else if (ioppbc .eq. 6) then
c       Truncated octahedron, axes normal to the square faces
        sum=dc1**2+dc2**2+dc3**2
        icell=1
        ad1=abs(dc1)
        ad2=abs(dc2)
        ad3=abs(dc3)
        dsum=ad1+ad2+ad3
        dmax=amax1(ad1,ad2,ad3)
        d1=-dmax*et1+et2
        d2=-dsum*et3+et4
        dminn=amin1(d1,d2)
        if (dminn .lt. 0.0) then
          sum=sum+dminn
          icell=2
        end if
      else if (ioppbc .eq. 6 .or. ioppbc .eq. 7) then
c       Truncated octahedron, X axis normal to a hexagon face
        icell=1
        d(1)=dc1
        d(2)=dc2
        d(3)=dc3
        call genimdist(d,cell,1,ncell,icmin,rmin2)
        if (icmin .gt. 1) icell=2
      else if (ioppbc .eq. 8) then
c       Hexagonal close-packed
        dax=abs(dc1)
        day=abs(dc2)
        daz=abs(dc3)
        dsx=dc1*dc1
        dsy=dc2*dc2
        dsz=dc3*dc3
        summ(1)=dsx+dsy+dsz
c       Distances from hexagonal plane neighbours
        dax2=(dax-edghalf)**2
        dax3=(dax-edge(1))**2
        day2=(day-edgsq3p2)**2
        summ(2)=dax2+day2+dsz
        summ(3)=dax3+dsy+dsz
c       Distances from face-touching upper neighbours
        day4=(dc2-edgpsq3)**2
        day5=(dc2+edgp2sq3)**2
        daz1=(daz-edgsq2p3)**2
        summ(4)=dsx+day4+daz1
        summ(5)=dax2+day5+daz1
c       Distances from vertex-touching upper neighbours
        day6=(dc2+edg2psq3)**2
        day7=(dc2-edgpsq3)**2
        summ(6)=dsx+day6+daz1
        summ(7)=dax3+day7+daz1
        sum=summ(1)
        icell=1
        do ic=2,7
          if (summ(ic) .lt. sum) then
            sum=summ(ic)
            icell=ic
           end if
        end do
      else if (ioppbc .eq. 0) then
c       Input pbc
        d(1)=dc1
        d(2)=dc2
        d(3)=dc3
        call genimdist(d,cell,1,ncell,icmin,sum)
      end if
c       Check
c       d(1)=dc1
c       d(2)=dc2
c       d(3)=dc3
c       call genimdist(d,cell,1,ncell,icmin,sumx)
c       if (abs(sum-sumx) .gt. 0.1 .or.
c    -           icmin .eq. 1 .and. icell  .gt. 1 .or.
c    -           icmin .gt. 1 .and. icell  .eq. 1)
c    -     print *,'DMC ERROR: sum,x=',sum,sumx,
c    -    ' icmin,icell=',icmin,icell
      return
      end
