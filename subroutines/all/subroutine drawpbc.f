      subroutine drawpbc(edgexyz,edge_gen,ioppbc,ixyzhex,cent,i2nd,
     -  icrot,crot,sizefac,ioutpdb)
      dimension edgexyz(3),edge_gen(3,3),cent(3),ixyzhex(3),crot(3,3),
     -  e2(3),vp(3,48),ixdup(48),hxy(2,6),vphl(3,7),vphu(3,7)
      common /graphics/ npixx,npixy,maxpixx,maxpixy,idwmain,idwplot,
     -  wx,wy,wz,wxdr
      common /rotmat/ matrot0(4,4),matrot(4,4),nomat0
      common /pbcrotmat/ torot_ac(3,3),torot_ca(3,3),tofac_ac,tofac_ca
      common /depthcuedat/ near,ifar,ramp0,idepth,idepthon,idrawh,
     -  linew,isidec,nbackb,idrawslv
c     print *,'DRAWPBC ioppbc=',ioppbc,' edge=',edgexyz,
c    -  ' sizefac=',sizefac
c     print *,'DRAWPBC cent=',cent
      if (ioppbc .gt. 9) return
      call indexit(ixdup,1,48,0)
      if (ioutpdb .gt. 0) iconntyp=2
      wxold=wx
      wx=0.0
c     Establish world coordinate system
      npixmin=min0(npixx,npixy)
      if (ioppbc .eq. 1) then
        edgemax=edgexyz(1)
        wx=(edgemax/2.0)*sizefac
        edgexyz(2)=edgexyz(1)
        edgexyz(3)=edgexyz(1)
      else if (ioppbc .eq. 2) then
c       Rectangular cell
        edgemax=amax1(edgexyz(1),edgexyz(2),edgexyz(3))
        wx=(edgemax/2.0)*sizefac
      else if (ioppbc .eq. 3) then
        wx=edgexyz(1)*sizefac
      else if (ioppbc .eq. 4) then
        wx=amax1(edgexyz(1)/2,edgexyz(2))*sizefac
      else if (ioppbc .eq. 5) then
        wx=amax1(edgexyz(1)/2,edgexyz(2)/2,edgexyz(3)/2)*sizefac
      else if (ioppbc .eq. 6 .or. ioppbc .eq. 7) then
        wx=edgexyz(1)*sizefac
      else if (ioppbc .eq. 8) then
        wx=(edgexyz(1)/sqrt(3.0))*sizefac
      else if (ioppbc .eq. 9) then
       wx=edge_gen(1,1)*1.3*sizefac
      else if (ioppbc .eq. 10) then
        print *,'Sorry, no cell is to be drawn for inputted images'
      else if (ioppbc .eq. 11) then
        print *,'Sorry, no cell is to be drawn for sphere PBC'
      else
        print *,'ERROR: invalid pbc code=',ioppbc
      end if
      wx=wx+amax1(abs(cent(1)),amax1(abs(cent(2)),abs(cent(3))))
      if (wx .lt. 5.0) wx=5.0
      wy=wx
      wz=wx
      if (wx .gt. wxold .and. i2nd .eq. 0) then
      end if
c     Start polyhedron object definition
      if (i2nd .eq. 0) then
c       Draw cell name
        wylab=0.9*wxdr
        wxlab=0.9*wxdr
c       print *,'wx,wxdr=',wx,wxdr,' wxlab,wylab=',wxlab,wylab
      end if
      if (ioppbc .eq. 1 .or. ioppbc .eq. 2) then
c       Rectangular cell
        nv=8
c       e2: half edges in pixel
        do k=1,3
          e2(k)=edgexyz(k)/2.0
        end do
c       vp(k,1-4): y-z face at the -x side
        vp(2,1)=-e2(2)
        vp(3,1)=-e2(3)
        vp(2,2)=+e2(2)
        vp(3,2)=-e2(3)
        vp(2,3)=+e2(2)
        vp(3,3)=+e2(3)
        vp(2,4)=-e2(2)
        vp(3,4)=+e2(3)
c       vp(k,5-8): y-z face at the +x side
        do m=1,4
          vp(1,m)=-e2(1)
          vp(1,m+4)=+e2(1)
          do k=2,3
            vp(k,m+4)=vp(k,m)
          end do
        end do
        call shiftmol(vp,8,cent,vp,1.0)
        call writevertices(vp,8,ioppbc,ioutpdb,iconntyp,icrot,crot)
        do m=1,4
c         Top polygon edge
          call connlinix(vp,m,mod(m,4)+1,ixdup,iconntyp,ioutpdb,48)
c         Bottom polygon edge
          call connlinix(vp,m+4,mod(m,4)+5,ixdup,iconntyp,ioutpdb,48)
c         Edge parallel to the prizm's axis
          call connlinix(vp,m,m+4,ixdup,iconntyp,ioutpdb,48)
        end do
      else if (ioppbc .eq. 3) then
c       FCC
        nv=14
        edge=edgexyz(1)
        edge2=edgexyz(1)/2.0
        call zeroit(vp,3*6)
        do k=1,3
          vp(k,2*(k-1)+1)=-edge
          vp(k,2*k)=edge
        end do
        do i=1,4
          vp(1,2*(i-1)+7)=(-1)**i*edge2
          vp(1,2*(i-1)+8)=vp(1,2*(i-1)+7)
          vp(2,2*(i-1)+7)=(-1)**((i-1)/2)*edge2
          vp(2,2*(i-1)+8)=vp(2,2*(i-1)+7)
          vp(3,2*(i-1)+7)=-edge2
          vp(3,2*(i-1)+8)=+edge2
        end do
        call shiftmol(vp,14,cent,vp,1.0)
        call finddup(vp,ixdup,14,nunique)
        call writevertices(vp,nunique,ioppbc,ioutpdb,iconntyp,
     -    icrot,crot)
        call connlinix(vp,1,07,ixdup,iconntyp,ioutpdb,48)
        call connlinix(vp,1,08,ixdup,iconntyp,ioutpdb,48)
        call connlinix(vp,1,11,ixdup,iconntyp,ioutpdb,48)
        call connlinix(vp,1,12,ixdup,iconntyp,ioutpdb,48)
        call connlinix(vp,2,09,ixdup,iconntyp,ioutpdb,48)
        call connlinix(vp,2,10,ixdup,iconntyp,ioutpdb,48)
        call connlinix(vp,2,13,ixdup,iconntyp,ioutpdb,48)
        call connlinix(vp,2,14,ixdup,iconntyp,ioutpdb,48)
        call connlinix(vp,3,11,ixdup,iconntyp,ioutpdb,48)
        call connlinix(vp,3,12,ixdup,iconntyp,ioutpdb,48)
        call connlinix(vp,3,13,ixdup,iconntyp,ioutpdb,48)
        call connlinix(vp,3,14,ixdup,iconntyp,ioutpdb,48)
        call connlinix(vp,4,07,ixdup,iconntyp,ioutpdb,48)
        call connlinix(vp,4,08,ixdup,iconntyp,ioutpdb,48)
        call connlinix(vp,4,09,ixdup,iconntyp,ioutpdb,48)
        call connlinix(vp,4,10,ixdup,iconntyp,ioutpdb,48)
        call connlinix(vp,5,07,ixdup,iconntyp,ioutpdb,48)
        call connlinix(vp,5,09,ixdup,iconntyp,ioutpdb,48)
        call connlinix(vp,5,13,ixdup,iconntyp,ioutpdb,48)
        call connlinix(vp,5,11,ixdup,iconntyp,ioutpdb,48)
        call connlinix(vp,6,08,ixdup,iconntyp,ioutpdb,48)
        call connlinix(vp,6,10,ixdup,iconntyp,ioutpdb,48)
        call connlinix(vp,6,14,ixdup,iconntyp,ioutpdb,48)
        call connlinix(vp,6,12,ixdup,iconntyp,ioutpdb,48)
      else if (ioppbc .eq. 4 .or. ioppbc .eq. 5) then
c       Hexagonal prism
        nv=12
        edgep2=edgexyz(1)/2.0
        if (ioppbc .eq. 4) then
          edgey=edgexyz(2)
          esq3p2=edgey*(sqrt(3.0)/2.0)
          vp(ixyzhex(2),1)=edgey
          vp(ixyzhex(2),2)=edgey/2.0
          vp(ixyzhex(3),2)=esq3p2
        else
          edgex=edgexyz(3)
          edgey=edgexyz(2)
          w=sqrt(edgey**2-edgex**2/4.0)
          h=(w-edgex/(2.0*sqrt(3.0)))/2.0
          vp(ixyzhex(2),1)=w-h
          vp(ixyzhex(3),2)=edgex/2.0
          vp(ixyzhex(2),2)=h
        end if
        vp(ixyzhex(3),1)=0.0
        vp(ixyzhex(3),3)=vp(ixyzhex(3),2)
        vp(ixyzhex(2),3)=-vp(ixyzhex(2),2)
        vp(ixyzhex(2),4)=-vp(ixyzhex(2),1)
        vp(ixyzhex(3),4)=0.0
        vp(ixyzhex(2),5)=-vp(ixyzhex(2),2)
        vp(ixyzhex(3),5)=-vp(ixyzhex(3),2)
        vp(ixyzhex(2),6)=-vp(ixyzhex(2),3)
        vp(ixyzhex(3),6)=-vp(ixyzhex(3),3)
        do m=1,6
          vp(ixyzhex(1),m)=-edgep2
          vp(ixyzhex(1),m+6)=+edgep2
          do k=2,3
            vp(ixyzhex(k),m+6)=vp(ixyzhex(k),m)
          end do
        end do
        call shiftmol(vp,12,cent,vp,1.0)
        call writevertices(vp,12,ioppbc,ioutpdb,iconntyp,icrot,crot)
        do m=1,6
c         Top polygon edge
          call connlinix(vp,m,mod(m,6)+1,ixdup,iconntyp,ioutpdb,48)
c         Bottom polygon edge
          call connlinix(vp,m+6,mod(m,6)+7,ixdup,iconntyp,ioutpdb,48)
c         Edge parallel to the prizm's axis
          call connlinix(vp,m,m+6,ixdup,iconntyp,ioutpdb,48)
        end do
      else if (ioppbc .eq. 6 .or. ioppbc .eq. 7) then
c       Truncated octahedron cell
        eto2=edgexyz(1)
        eto=edgexyz(1)/2.0
        print *,'eto,sto2=',eto,eto2
        call zeroit(vp,(3*6*8))
c       Hexagon face (+x,+y,+z) quadrant
        vp(1,1)=eto
        vp(3,1)=eto2
        vp(2,2)=eto
        vp(3,2)=eto2
        vp(2,3)=eto2
        vp(3,3)=eto
        vp(1,4)=eto
        vp(2,4)=eto2
        vp(1,5)=eto2
        vp(2,5)=eto
        vp(1,6)=eto2
        vp(3,6)=eto
c       Hexagon face 2 (+x,+y,-z) quadrant
        vp(1,1+6)=eto
        vp(3,1+6)=-eto2
        vp(2,2+6)=eto
        vp(3,2+6)=-eto2
        vp(2,3+6)=eto2
        vp(3,3+6)=-eto
        vp(1,4+6)=eto
        vp(2,4+6)=eto2
        vp(1,5+6)=eto2
        vp(2,5+6)=eto
        vp(1,6+6)=eto2
        vp(3,6+6)=-eto
c       Hexagon face 3 (-x,+y,-z) quadrant
        vp(1,1+12)=-eto
        vp(3,1+12)=-eto2
        vp(2,2+12)=eto
        vp(3,2+12)=-eto2
        vp(2,3+12)=eto2
        vp(3,3+12)=-eto
        vp(1,4+12)=-eto
        vp(2,4+12)=eto2
        vp(1,5+12)=-eto2
        vp(2,5+12)=eto
        vp(1,6+12)=-eto2
        vp(3,6+12)=-eto
c       Hexagon face 4 (-x,+y,+z) quadrant
        vp(1,1+18)=-eto
        vp(3,1+18)=eto2
        vp(2,2+18)=eto
        vp(3,2+18)=eto2
        vp(2,3+18)=eto2
        vp(3,3+18)=eto
        vp(1,4+18)=-eto
        vp(2,4+18)=eto2
        vp(1,5+18)=-eto2
        vp(2,5+18)=eto
        vp(1,6+18)=-eto2
        vp(3,6+18)=eto
c       Hexagon face 5 (+x,-y,+z) quadrant
        vp(1,1+24)=eto
        vp(3,1+24)=eto2
        vp(2,2+24)=-eto
        vp(3,2+24)=eto2
        vp(2,3+24)=-eto2
        vp(3,3+24)=eto
        vp(1,4+24)=eto
        vp(2,4+24)=-eto2
        vp(1,5+24)=eto2
        vp(2,5+24)=-eto
        vp(1,6+24)=eto2
        vp(3,6+24)=eto
c       Hexagon face 6 (+x,-y,-z) quadrant
        vp(1,1+30)=eto
        vp(3,1+30)=-eto2
        vp(2,2+30)=-eto
        vp(3,2+30)=-eto2
        vp(2,3+30)=-eto2
        vp(3,3+30)=-eto
        vp(1,4+30)=eto
        vp(2,4+30)=-eto2
        vp(1,5+30)=eto2
        vp(2,5+30)=-eto
        vp(1,6+30)=eto2
        vp(3,6+30)=-eto
c       Hexagon face 7 (-x,-y,-z) quadrant
        vp(1,1+36)=-eto
        vp(3,1+36)=-eto2
        vp(2,2+36)=-eto
        vp(3,2+36)=-eto2
        vp(2,3+36)=-eto2
        vp(3,3+36)=-eto
        vp(1,4+36)=-eto
        vp(2,4+36)=-eto2
        vp(1,5+36)=-eto2
        vp(2,5+36)=-eto
        vp(1,6+36)=-eto2
        vp(3,6+36)=-eto
c       Hexagon face 8 (-x,-y,+z) quadrant
        vp(1,1+42)=-eto
        vp(3,1+42)=eto2
        vp(2,2+42)=-eto
        vp(3,2+42)=eto2
        vp(2,3+42)=-eto2
        vp(3,3+42)=eto
        vp(1,4+42)=-eto
        vp(2,4+42)=-eto2
        vp(1,5+42)=-eto2
        vp(2,5+42)=-eto
        vp(1,6+42)=-eto2
        vp(3,6+42)=eto
        if (ioppbc .eq. 7)
     -    call rotate_c(vp,48,torot_ca,vp,'VERTICES',8)
        call shiftmol(vp,48,cent,vp,1.0)
        call finddup(vp,ixdup,48,nunique)
        call writevertices(vp,nunique,ioppbc,ioutpdb,iconntyp,
     -    icrot,crot)
        do ifa=1,8
          do m=1,6
            call connlinix(vp,m+(ifa-1)*6,mod(m,6)+1+(ifa-1)*6,ixdup,
     -        iconntyp,ioutpdb,48)
          end do
        end do
      else if (ioppbc .eq. 8) then
C       Hexagonal close-packed
        d=edgexyz(1)
        t=d*sqrt(3.0)/(2.0*sqrt(2.0))
        hl=d/(2.0*sqrt(6.0))
        hh=(t+hl)/2.0
c       print *,'t,hh,hl=',t,hh,hl
c       x-y coordinates of wrapping hexagon vertices
        hxy(1,1)=d/2.0
        hxy(2,1)=-d/(2.0*sqrt(3.0))
        hxy(1,2)=0.0
        hxy(2,2)=-d/sqrt(3.0)
        hxy(1,3)=-d/2.0
        hxy(2,3)=-d/(2.0*sqrt(3.0))
        hxy(1,4)=-d/2.0
        hxy(2,4)=+d/(2.0*sqrt(3.0))
        hxy(1,5)=0.0
        hxy(2,5)=d/sqrt(3.0)
        hxy(1,6)=d/2.0
        hxy(2,6)=+d/(2.0*sqrt(3.0))
        call zeroit(vp,6)
        vp(3,1)=t
        vp(3,2)=-t
c       Top rhombuses
c       First the chair hexagon
        do i=1,6
          call trnsfr(vphu(1,i),hxy(1,i),2)
          if (mod(i,2) .eq. 1) then
            vphu(3,i)=hl
          else
            vphu(3,i)=hh
          end if
        end do
        call trnsfr(vp(1,3),vphu,18)
c       Bottom rhombuses
        call trnsfr(vphl(1,1),vphu(1,1),21)
        do k=1,7
          vphl(3,k)=-vphl(3,k)
        end do
        call trnsfr(vp(1,9),vphl,18)
        call finddup(vp,ixdup,14,nunique)
        call writevertices(vp,nunique,ioppbc,ioutpdb,iconntyp,
     -    icrot,crot)
c       Top rhombuses
        do m=1,6
          call connlinix(vp,m+2,mod(m,6)+3,ixdup,iconntyp,ioutpdb,48)
        end do
        do k=2,6,2
          call connlinix(vp,1,k+2,ixdup,iconntyp,ioutpdb,48)
        end do
c       Bottom rhombuses
        do m=1,6
          call connlinix(vp,8+m,mod(m,6)+9,ixdup,iconntyp,ioutpdb,48)
        end do
        do k=2,6,2
          call connlinix(vp,2,k+8,ixdup,iconntyp,ioutpdb,48)
        end do
c       Side connections
        do i=1,6
          call connlinix(vp,i+2,i+8,ixdup,iconntyp,ioutpdb,48)
        end do
      else if (ioppbc .eq. 9) then
c       Octahedral - draw parallelepiped
c       print *,'DRAWPBC ioppbc=9 '
        do k=1,3
          vp(k,1)=(-edge_gen(k,1)-edge_gen(k,2)-edge_gen(k,3))/2.0
        end do
        call arrsum(vp(1,1),edge_gen(1,1),vp(1,2),3)
        call arrsum(vp(1,2),edge_gen(1,2),vp(1,3),3)
        call arrsum(vp(1,1),edge_gen(1,2),vp(1,4),3)
        call arrsum(vp(1,1),edge_gen(1,3),vp(1,5),3)
        call arrsum(vp(1,5),edge_gen(1,1),vp(1,6),3)
        call arrsum(vp(1,6),edge_gen(1,2),vp(1,7),3)
        call arrsum(vp(1,5),edge_gen(1,2),vp(1,8),3)
        call writevertices(vp,8,ioppbc,ioutpdb,iconntyp,
     -    icrot,crot)
        do m=1,4
          call connlinix(vp,m,mod(m,4)+1,ixdup,iconntyp,ioutpdb,48)
          call connlinix(vp,m+4,mod(m,4)+5,ixdup,iconntyp,ioutpdb,48)
          call connlinix(vp,m,m+4,ixdup,iconntyp,ioutpdb,48)
        end do
      end if
c        write (6,7171) (i,(vp(k,i),k=1,3),i=1,nv)
c7171    format(i5,3f10.3)
      return
      end
