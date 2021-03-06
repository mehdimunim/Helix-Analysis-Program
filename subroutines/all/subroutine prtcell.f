      subroutine prtcell(ioppbc,edge,edge_gen,r,vol,nw,iprtopt)
      dimension edge(3),edge_gen(3,3)
      common /numbers/ sq3,sq3inv,sq3p2,sq2p3
      dimension c0(3)
      if (iprtopt .lt. 0) then
        iprt=3
        iout=-iprtopt
      else
        iprt=iprtopt
        iout=6
      end if
      if (ioppbc .eq. 1) then
        vol=edge(1)**3
        if (iprt .gt. 0) write (iout,1001) edge(1)
      else if (ioppbc .eq. 2) then
        vol=edge(1)*edge(2)*edge(3)
        if (iprt .gt. 0) write (iout,1002) edge
      else if (ioppbc .eq. 3) then
        vol=2.0*edge(1)**3
        if (iprt .gt. 0) write (iout,1003) edge(1)
      else if (ioppbc .eq. 4) then
        vol=edge(1)*edge(2)**2*3.0*sqrt(3.0)/2.0
        if (iprt .gt. 0) write (iout,1004) edge(1),edge(2)
      else if (ioppbc .eq. 5) then
        edgex=edge(3)
        edgey=edge(2)
        edgep=edge(1)
        w=sqrt(edgey**2-edgex**2/4.0)
        h=(w-edgex/(2.0*sqrt(3.0)))/2.0
        vol=edgep*edgex*w
        obang=180.0-(180.0/3.14159)*acos(edgex/(2.0*edgey))
        if (iprt .gt. 0) write (iout,1005) edgep,obang,edgey,edgex
      else if (ioppbc .eq. 6 .or. ioppbc .eq. 7) then
        vol=4.0*edge(1)**3
        if (iprt .gt. 0) then
          rhex=edge(1)/sq3p2
          cx_c=2.0*edge(1)
          cx_a=cx_c*sq3p2
          write (iout,1006) edge(1),rhex,cx_c,cx_a
        end if
      else if (ioppbc .eq. 8) then
        vol=edge(1)**3/sqrt(2.0)
        if (iprt .gt. 0) write (iout,1008) edge(1)
      else if (ioppbc .eq. 9) then
        call vprd(edge_gen(1,2),edge_gen(1,3),c0)
        vol=abs(scprod(edge_gen(1,1),c0))
        if (iprt .gt. 0) write (iout,1009) edge(1)
      else if (ioppbc .eq. 10) then
        write (iout,*) 'Inputtted image cell centers'
      else if (ioppbc .eq. 11) then
        vol=4.0*3.141592/3.0*r**3
        iprt=3
        if (iprt .gt. 0) write (iout,1010) r
      else if (ioppbc .eq. -1) then
        vol=4.0*3.141592/3*r**3
        iprt=3
      else
        print *,'PROGRAM ERROR: invalid PBC option in prtcell:',ioppbc
        return
      end if
      nw=vol/(30.090)
      if (iprt .eq. 1) write (iout,2000) vol,nw
      if (iprt .eq. 2) write (iout,2001) vol,nw
      return
1001  format(' Edge of the cube=',f11.6,' A')
1002  format(' Edges of the rectangle=',3f11.6,' A')
1003  format(' FCC cell parameter used=',f11.6,' A')
1004  format(' Hexagonal prism length=',f11.6,
     -  ' A hexagon edge=',f11.6,' A')
1005  format(' Hexagonal prism (skewed) length=',f11.6,
     -  ' A Obtuse angle=',f8.2,' deg',/,' Obtuse axis length=',f11.6,
     -  ' A Cartesian axis length=',f11.6,' A')
1006  format(' Truncated octahedron',/,
     -  ' Distance of square face from the center=',f11.6,' A',/,
     -  ' Distance of hexagonal face from the center=',f11.6,' A',/,
     -  ' Charmm     cell X-parameter=',f11.6,' A',/,
     -  ' Amber/NAMD cell X-parameter=',f11.6,' A')
1008  format(' Diameter of the close-packed sphere=',f11.6,' A')
1009  format(' Edge length of the paralellepid=',f11.6,' A')
1010  format(' Radius of the boundary sphere=',f11.6,' A')
2000  format(' Volume=',f15.2,' A**3',' Number of waters that fit in=',
     -  i9)
2001  format(' Optimal volume=',f11.2,' A**3',
     -  ' Number of waters that fit in=',i6)
      end
