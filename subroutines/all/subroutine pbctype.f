      subroutine pbctype(ioppbc,npbc,ixyzhex,nonone)
      dimension ixyzhex(3)
      common /pbcrotmat/ torot_ac(3,3),torot_ca(3,3),tofac_ac,tofac_ca
      common /rotmat/ matrot0(4,4),matrot(4,4),nomat0
      character*1 xyz
      common /axislab/ xyz(3)
      character*1 ans,ansedge,ansvertex
100   call quiz(ans,ians,'r',' ',0,'PBC cell type',13,0,5,6,0)
      if (ans .eq. 'n') then
c       None
        if (nonone .eq. 0) then
          ioppbc=-1
          npbc=1
        else
          print *,'No PBC choice is invalid in this context'
          go to 100
        end if
      else if (ans .eq. 'c') then
c       Cubic
        ioppbc=1
        npbc=1
      else if (ans .eq. 'r') then
c       Rectangular
        ioppbc=2
        npbc=3
      else if (ans .eq. 'f') then
c       Face-centered cubic
        ioppbc=3
        npbc=1
c       Face-centered cubic
      else if (ans .eq. 'h' .or. ans .eq. 'w') then
c       Hexagonal prizms
200     call quiz(ansedge,ixyzhex(1),' ',
     -    'parallel to the prism length',-28,'axis',4,0,5,6,0)
        call quiz(ansvertex,ixyzhex(2),' ',
     -    'going through the hexagon vertex',-32,'axis',4,0,5,6,0)
        if (ansedge .eq. ansvertex) then
          print *,'ERROR: can not use the same axes twice'
          go to 200
        end if
        ixyzhex(3)=1
        do while (ixyzhex(3) .eq. ixyzhex(1) .or.
     -    ixyzhex(3) .eq. ixyzhex(2))
          ixyzhex(3)=ixyzhex(3)+1
        end do
        write (6,2000) xyz(ixyzhex(1)),xyz(ixyzhex(2))
        if (ans .eq. 'h') then
          ioppbc=4
          npbc=2
        else
c         Skewed hexagonal prizm
          ioppbc=5
          npbc=3
        end if
      else if (ans .eq. 't' .or. ans .eq. 'o') then
c       Truncated Octahedron,
        npbc=1
        torot_ac(1,1)=1.0/sqrt(3.0)
        torot_ac(2,1)=1.0/sqrt(3.0)
        torot_ac(3,1)=1.0/sqrt(3.0)
        torot_ac(1,2)=-1.0/sqrt(6.0)
        torot_ac(2,2)=-1.0/sqrt(6.0)
        torot_ac(3,2)=2.0/sqrt(6.0)
        torot_ac(1,3)=1.0/sqrt(2.0)
        torot_ac(2,3)=-1.0/sqrt(2.0)
        torot_ac(3,3)=0.0
        tofac_ac=2.0/sqrt(3.0)
        torot_ca(1,1)=1.0/sqrt(3.0)
        torot_ca(1,2)=1.0/sqrt(3.0)
        torot_ca(1,3)=1.0/sqrt(3.0)
        torot_ca(2,1)=-1.0/sqrt(6.0)
        torot_ca(2,2)=-1.0/sqrt(6.0)
        torot_ca(2,3)=2.0/sqrt(6.0)
        torot_ca(3,1)=1.0/sqrt(2.0)
        torot_ca(3,2)=-1.0/sqrt(2.0)
        torot_ca(3,3)=0.0
        tofac_ca=sqrt(3.0)/2.0
        if (ans .eq. 't') then
c         Charmm convention (x normal to square)
          ioppbc=6
          write (6,2011)
        else if (ans .eq. 'o') then
c         Truncated Octahedron, Amber/NAMD convention (X normal to hexagon)
          ioppbc=7
          npbc=1
          write (6,2012)
        end if
      else if (ans .eq. 'x') then
c       Hexagonal close packed
        ioppbc=8
        npbc=1
      else if (ans .eq. 'd') then
c       Octahedral
        ioppbc=9
        npbc=1
        write (6,2013)
      else if (ans .eq. 'i') then
c       Image-cell file input
        ioppbc=10
      else if (ans .eq. 's') then
c       Spehere
        ioppbc=11
        npbc=1
      else
        print *,'PROGRAM ERROR'
      end if
      return
2000  format(' Hexagonal prism edge: along the ',a,' axis',/,
     -  ' Axis going through the hexagon: ',a)
2011  format(' Charmm convention - X axis is normal to a square')
2012  format(' Amber/NAMD convention - X axis is normal to a hexagon')
2013  format(' NOTE: For now, this PBC can not be used for optimizing',
     -  ' orientation')
      end
