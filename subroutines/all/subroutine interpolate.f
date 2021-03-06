      subroutine interpolate(x,y,z,gx,gy,gz,xstart,ystart,zstart,phi)
      parameter (MAXPHI=400)
      common /nnwork/ phimap(MAXPHI,MAXPHI,MAXPHI)
      ix=(x-xstart)/gx+1
      iy=(y-ystart)/gy+1
      iz=(z-zstart)/gz+1
      xm=xstart+(ix-1)*gx
      ym=ystart+(iy-1)*gy
      zm=zstart+(iz-1)*gz
      wxp=(x-xm)/gx
      wyp=(y-ym)/gy
      wzp=(z-zm)/gz
      wxm=1.0-wxp
      wym=1.0-wyp
      wzm=1.0-wzp
      p0y0=wym*phimap(ix+0,iy+0,iz+0)+wyp*phimap(ix+0,iy+1,iz+0)
      p0y1=wym*phimap(ix+0,iy+0,iz+1)+wyp*phimap(ix+0,iy+1,iz+1)
      p1y0=wym*phimap(ix+1,iy+0,iz+0)+wyp*phimap(ix+1,iy+1,iz+0)
      p1y1=wym*phimap(ix+1,iy+0,iz+1)+wyp*phimap(ix+1,iy+1,iz+1)
      p0yz=wzm*p0y0+wzp*p0y1
      p1yz=wzm*p1y0+wzp*p1y1
      phi=wxm*p0yz+wxp*p1yz
      return
      end
