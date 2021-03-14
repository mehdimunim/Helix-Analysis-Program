      subroutine crorgn(edg,edge_gen,ioppbc,ndim,ncell,cell,cellalt,
     -  ixyzhex,rinscr,rcirc)
c#    MMC routine 064 lstmod: 09/12/90
c*****Generate the negative cell-center coordinates for the pbc routines
      dimension edg(3),edge_gen(3,3),cell(3,27),cellalt(3,27)
      common /pbcrotmat/ torot_ac(3,3),torot_ca(3,3),tofac_ac,tofac_ca
      dimension y(3),x(3),yx(3),yy(3),ixyzhex(3),cx(3),cy(3),cz(3)
      data edgex /0.0/,edgey /0.0/,edge /0.0/,edgep /0.0/,d /0.0/
      if (ioppbc .eq. -1) return
c     print *,'CRORGN ioppbc=',ioppbc,' edg=',edg
      onethird=1.0/3.0
      sq22p3=sqrt(2.0)*2.0/3.0
      sq3=sqrt(3.0)
      sq3p2=sq3/2.0
      sq2p3=sqrt(2.0/3.0)
      if (ioppbc .eq. 1) then
c       Cubic
        edgex=edg(1)
        edgey=edg(1)
        edge=edg(1)
      else if (ioppbc .eq. 2) then
c       Rectangular
        edgex=edg(1)
        edgey=edg(2)
        edge=edg(3)
      else if (ioppbc .eq. 3) then
c       Face-centered cubic
        edge=edg(1)
      else if (ioppbc .eq. 4) then
c       Hexagonal prizm (regular)
        edgex=edg(1)
        edge=edg(2)
      else if (ioppbc .eq. 5) then
c       Hexagonal prizm (skewed)
        edgep=edg(1)
        edgex=edg(3)
        edgey=edg(2)
      else if (ioppbc .eq. 6 .or. ioppbc .eq. 7) then
c       Truncated octahedron
cDAN    edge=edg(1)*2.0*sqrt(2.0)
        edge=edg(1)
      else if (ioppbc .eq. 8) then
c       Hexagonal close packed
        d=edg(1)
      else if (ioppbc .eq. 9) then
c       Octahedral
        call zeroit(cx,3)
        cx(1)=edg(1)
c       do k=1,3
c         cy(k)=cx(k)*onethird-edge_gen(k,2)*sq22p3
c       end do
c       solve (ez.ex)=e**2/3; (ez.ey)=-e**2/3; (ez.ez)=e**2
c       cz(3)=edg(1)
c       cz(1)=(edg(1)**2/3.0*(cy(2)+cx(2))+edg(1)*(cy(3)-cx(3)))/
c    -    (cx(1)*cy(2)-cy(1)*cx(2))
c       cz(2)=(-edg(1)**2/3.0-cz(1)*cy(1)-cz(3)*cy(3))/cy(2)
c       csum=cz(1)**2+cz(2)**2+cz(3)**2
c       fac=edg(1)/sqrt(csum)
c       do k=1,3
c         cz(k)=cz(k)*fac
c       end do
        cy(1)=-edg(1)/3.0
        cy(2)=edg(1)*sq22p3
        cy(3)=0.0
        cz(1)=edg(1)/3.0
        cz(2)=cy(2)/2.0
        cz(3)=edg(1)*sq2p3
        call trnsfr(edge_gen(1,1),cx,3)
        call trnsfr(edge_gen(1,2),cy,3)
        call trnsfr(edge_gen(1,3),cz,3)
c        write (6,9877) edge_gen
      else if (ioppbc .eq. 11) then
c       Sphere
        d=edg(1)
      end if
      threp2=1.5
c     cell will contain the negative of the respective cell center coords
      call zeroit(cell,3*27)
c     print *,'IOPPBC=',ioppbc
      go to (10,10,30,40,50,60,60,80,90,100,110),ioppbc
      print *,'PROGRAM ERROR: invalid ioppbc in crorgn=',ioppbc
      stop
c-----Siple cubic or rectangular
10    ncell=27
      y(1)=0.
      y(2)=edge
      y(3)=- edge
      yx(1)=0.0
      yx(2)=edgex
      yx(3)=-edgex
      yy(1)=0.0
      yy(2)=edgey
      yy(3)=-edgey
      ii=0
      do i=1,ndim
        x(3)=y(i)
        do k=1,ndim
          x(2)=yy(k)
          do l=1,ndim
            x(1)=yx(l)
            ii=ii+1
            do j=1,3
              cell(j,ii)=x(j)
            end do
          end do
        end do
      end do
      rcirc=sqrt(edgex**2+edgey**2+edge**2)/2.0
      rinscr=amin1(edgex,edgey,edge)/2
      go to 999
c-----Face-centered cubic
30    ncell=19
      do i=2,5
        cell(2,i)=-edge
        cell(3,i)=-edge
      end do
      cell(2,2)=edge
      cell(3,2)=edge
      cell(2,4)=edge
      cell(3,5)=edge
      do i=2,5
        cell(1,4+i)=cell(2,I)
        cell(3,4+i)=cell(3,i)
        cell(1,8+i)=cell(2,I)
        cell(2,8+i)=cell(3,i)
      end do
      cell(1,14)=-2.0*edge
      cell(1,15)=-cell(1,14)
      cell(2,16)=cell(1,14)
      cell(2,17)=-cell(2,16)
      cell(3,18)=cell(1,14)
      cell(3,19)=-cell(3,18)
      rcirc=edge
      rinscr=edge/sqrt(2.0)
      go to 999
c-----Hexagonal prism (regular)
40    continue
      ncell=21
c     Vertex of the hexagon along the ixyzhex(2) axis
      cell(ixyzhex(3),2)=-edge*sq3p2
      cell(ixyzhex(2),2)=-edge*threp2
      cell(ixyzhex(3),3)=-edge*sq3
      cell(ixyzhex(3),4)=-edge*sq3p2
      cell(ixyzhex(2),4)=edge*threp2
      do k=5,7
        l=9-k
        cell(ixyzhex(2),k)=cell(ixyzhex(2),l)
        cell(ixyzhex(3),k)=-cell(ixyzhex(3),l)
      end do
      do k=1,7
        l1=k+7
        l2=k+14
        cell(ixyzhex(1),l1)=edgex
        cell(ixyzhex(1),l2)=-edgex
        cell(ixyzhex(2),l1)=cell(ixyzhex(2),k)
        cell(ixyzhex(2),l2)=cell(ixyzhex(2),k)
        cell(ixyzhex(3),l1)=cell(ixyzhex(3),k)
        cell(ixyzhex(3),l2)=cell(ixyzhex(3),k)
      end do
      rcirc=sqrt((edgex/2.0)**2+edge**2)
      rinscr=amin1(edgex/2.0,edge*sq3p2)
      go to 999
c-----Hexagonal prism (skewed)
50    continue
      ncell=21
      w=sqrt(edgey**2-edgex**2/4.0)
      h=(w-edgex/(2.0*sqrt(3.0)))/2.0
c     Vertex of the hexagon along the ixyzhex(2) axis
      cell(ixyzhex(3),2)=-edgex/2.0
      cell(ixyzhex(2),2)=-w
      cell(ixyzhex(3),3)=-edgex
      cell(ixyzhex(3),4)=cell(ixyzhex(3),2)
      cell(ixyzhex(2),4)=-cell(ixyzhex(2),2)
      do k=5,7
        l=9-k
        cell(ixyzhex(2),k)=cell(ixyzhex(2),l)
        cell(ixyzhex(3),k)=-cell(ixyzhex(3),l)
      end do
      do k=1,7
        l1=k+7
        l2=k+14
        cell(ixyzhex(1),l1)=edgep
        cell(ixyzhex(1),l2)=-edgep
        cell(ixyzhex(2),l1)=cell(ixyzhex(2),k)
        cell(ixyzhex(2),l2)=cell(ixyzhex(2),k)
        cell(ixyzhex(3),l1)=cell(ixyzhex(3),k)
        cell(ixyzhex(3),l2)=cell(ixyzhex(3),k)
      end do
      rcirc=sqrt((edgep/2.0)**2+amax1((w-h)**2,(edgex/2.0)**2+h**2))
      rinscr=amin1(edgep,edgex,w-h)/2.0
      go to 999
c-----Truncated octahedon
60    ncell=15
c     Truncated face transforms
c     print *,'CRORG edge=',edge,' ncell=',ncell
      do i=2,6,2
        cell(i/2,i) =2.0*edge
        cell(i/2,i+1) = -2.0*edge
      end do
c     +/-z, xy face transforms
      ic=7
      do ix=1,2
        do iy=1,2
          do iz=1,2
            ic=ic+1
            cell(1,ic)=(-1)**(ix)*edge
            cell(2,ic)=(-1)**(iy)*edge
            cell(3,ic)=(-1)**(iz)*edge
          end do
        end do
      end do
      if (ioppbc .eq. 7) then
        call trnsfr(cellalt,cell,ncell*3)
        call rotate_c(cell,ncell,torot_ca,cell,'CRORGN',6)
      end if
      rcirc=edge*sqrt(5.0)/2.0
      rinscr=edge
      go to 999
c-----Hexagonal close packed
c     Neighbour cells in the hexagonal plane
80    ncell=19
      cell(1,2)=d
      cell(1,3)=d/2.0
      cell(2,3)=d*sq3p2
      cell(1,4)=d/2.0
      cell(2,4)=-d*sq3p2
      cell(1,5)=-d
      cell(1,6)=-d/2.0
      cell(2,6)=d*sq3p2
      cell(1,7)=-d/2.0
      cell(2,7)=-d*sq3p2
c     Neighbour cells above with touching face
      cell(2,8)=d/sq3
      cell(1,9)=d/2.0
      cell(2,9)=-d/(2.0*sq3)
      cell(1,10)=-d/2.0
      cell(2,10)=-d/(2.0*sq3)
c     Neighbour cells above with touching vertex
      cell(2,11)=-d*2.0/sq3
      cell(1,12)=d
      cell(2,12)=+d/sq3
      cell(1,13)=-d
      cell(2,13)=d/sq3
      call trnsfr(cell(1,14),cell(1,8),18)
      do k=8,13
        cell(1,k+6)=cell(1,k)
        cell(2,k+6)=cell(2,k)
        cell(3,k)=d*sq2p3
        cell(3,k+6)=-cell(3,k)
      end do
      rcirc=d*3.0/(sqrt(5.0)*2.0)
      rinscr=d/2.0
      go to 999
c-----Octahedral
90    ncell=1
c      write (6,9877) edge_gen
c9877  format(' E_G: ',3(3f8.4,1x))
      do ix=-1,1
        do iy=-1,1
          do iz=-1,1
            if (ix .ne. 0 .or. iy .ne. 0 .or. iz .ne. 0) then
              ncell=ncell+1
              do k=1,3
                cell(k,ncell)=ix*edge_gen(k,1)+iy*edge_gen(k,2)+
     -            iz*edge_gen(k,3)
              end do
            end if
          end do
        end do
      end do
      if (ncell .ne. 27) then
        print *,'PROGRAM ERROR: ncell=',ncell,' (instead of 27)'
        stop
      end if
      do k=1,3
        cz(k)=(edge_gen(k,1)+edge_gen(k,2)+edge_gen(k,3))/2.0
      end do
      rinscr=edg(1)*onethird/2.0
      rcirc=sqrt(scprod(cz,cz))
      go to 999
c-----Input cell centers
100   ncell=1
      rinscr=0.0
      rcirc=0.0
      go to 999
c-----Sphere
110   ncell=1
      rinscr=d
      rcirc=d
999   return
      end
