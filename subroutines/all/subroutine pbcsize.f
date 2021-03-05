      subroutine pbcsize(ioppbc,edge,npbc)
      dimension edge(3)
      character*1 xyz
      common /axislab/ xyz(3)
      common /numbers/ sq3,sq3inv,sq3p2,sq2p3
      if (ioppbc .eq. 0) then
c       Will be read with the images
      else if (ioppbc .eq. 1) then
        call getreal('Edge length (A)',15,999999.0,edge(1),1,0)
      else if (ioppbc .eq. 2) then
        do k=1,3
          call getreal('Edge length in the '//xyz(k)//' direction (A)',
     -      34,999999.0,edge(k),1,0)
        end do
      else if (ioppbc .eq. 3) then
        call getreal('Edge parameter of the FCC cell (A)',34,999999.0,
     -    edge(1),1,58)
      else if (ioppbc .eq. 4) then
c       Hexagonal prism (regular)
        call getreal('Length of the prism (A)',23,999999.0,edge(1),1,0)
        call getreal('Edge of the hexagon (A)',23,999999.0,edge(2),1,0)
      else if (ioppbc .eq. 5) then
c       Skewed hexagonal prism
        call getreal('Length of the prism (A)',23,999999.0,edge(1),1,0)
        call getreal('Cell length (a) along Cartesian axis (A)',40,
     -    999999.0,edge(3),1,59)
        call getreal(
     -    'Cell length (b) at ca 120 deg of Cartesian axis (A)',51,
     -    999999.0,edge(2),1,59)
      else if (ioppbc .eq. 6) then
c       Truncated octahedron, Charmm convention (x axis to a square face)
        call getreal('Charmm periodic cell X coordinate (A)',37,
     -    999999.0,cellx,1,60)
        edge(1)=cellx/2.0
      else if (ioppbc .eq. 7) then
c       Truncated octahedron, Amber/NAMD convention (x axis to a hexagon face)
        call getreal('Amber/NAMD periodic cell X coordinate (A)',41,
     -    999999.0,cellx,1,60)
        edge(1)=(cellx/2.0)/sq3p2
      else if (ioppbc .eq. 8) then
c       Hexagonal close packing
        call getreal('Inscribed sphere diameter (A)',29,999999.0,
     -    edge(1),1,61)
      else if (ioppbc .eq. 9) then
        call getreal('Edge of the rhomboid (A)',24,999999.0,
     -    edge(1),1,00)
      else if (ioppbc .eq. 10) then
        call getreal('Radius of the sphere',20,999999.0,edge(1),1,0)
      end if
      do k=npbc+1,3
        edge(k)=edge(npbc)
      end do
      return
      end
