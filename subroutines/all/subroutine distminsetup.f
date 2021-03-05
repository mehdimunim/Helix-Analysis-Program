      subroutine distminsetup(edge,ioppbc)
      dimension edge(3)
      common /distmindata/ ef1,ef2,ef3,ef4,edgep,edgep2,edgey,edgex,
     -  edgex2,w,h,et1,et2,et3,et4,edge2(3),edghalf,edgsq3p2,edgp2sq3,
     -  edgsq2p3,edgpsq3,edg2psq3
      common /numbers/ sq3,sq3inv,sq3p2,sq2p3
      if (ioppbc .eq. 3) then
c       Face-centered cubic
        ef1=2.0*edge(1)
        ef2=4.0*edge(1)
        ef3=2.0*edge(1)**2
        ef4=4.0*edge(1)**2
      else if (ioppbc .eq. 4) then
c       Regular hexagonal prizm
        edgep=edge(1)
        edgep2=edgep/2.0
        w=edge(2)*1.5
        h=edge(2)/2.0
        edgex2=edge(2)*sq3p2
        edgex=edge(2)*sq3
      else if (ioppbc .eq. 5) then
c       Skewed hexagonal prizm
        edgep=edge(1)
        edgep2=edgep/2.0
        edgey=edge(2)
        edgex=edge(3)
        edgex2=edgex/2.0
        w=sqrt(edgey**2-edgex**2/4.0)
        h=(w-edgex/(2.0*sq3))/2.0
      else if (ioppbc .eq. 6 .or. ioppbc .eq. 7) then
c       Truncated octahedron
        et1=4.0*edge(1)
        et2=4.0*edge(1)**2
        et3=2.0*edge(1)
        et4=3.0*edge(1)**2
      else if (ioppbc .eq. 1) then
c       Cubic
        do k=1,3
          edge(k)=edge(1)
          edge2(k)=edge(1)/2.0
        end do
      else if (ioppbc .eq. 8) then
c       Hexagonal close-packed
        edghalf=edge(1)/2.0
        edgsq3p2=edge(1)*sq3p2
        edgsq2p3=edge(1)*sqrt(2.0/3.0)
        edgp2sq3=edge(1)/(2.0*sq3)
        edgpsq3=edge(1)/sq3
        edg2psq3=edge(1)/sq3p2
      else
c       Rectangular
        do k=1,3
          edge2(k)=edge(k)/2.0
        end do
      end if
      return
      end
