      subroutine getedge(co,cn,ih,n,nnh,c0,edge,edge_gen,ioppbc,npbc,
     -  ireduce,cell,ncell,cellalt,ixyzhex,closlim,list,iatnum,cv,
     -  maxrec)
      dimension list(maxrec),iatnum(maxrec)
      dimension co(3,n),cn(3,n),ih(n),c0(3),edge(3),edge_gen(3,3),
     -  cell(3,27),cv(n),cellalt(3,27),ixyzhex(3)
      character*1 xyz
      common /axislab/ xyz(3)
      dimension cmin(3),cmax(3),edgen(3)
      data kmax /0/
c     Calculate the cell parameters from the molecule size and the solvent
c     layer width
      width=closlim/2.0
      if (ioppbc .eq. 3 .or. (ioppbc .ge. 6 .and. ioppbc .le. 8)) then
c       Mostly spherical PBC's (FCC, TO or HEX CP)
        print *,'Cell parameter obtained from smallest enclosing ',
     -          'sphere:'
        call compact(co,cn,iatnum,ih,n,nnh,c0,rmin,rorgext,rorgcom,
     -    list,np,nlist,1,cv,maxrec)
        write (6,1003) rmin
        if (ioppbc .eq. 3) then
c         FCC
          edge(1)=(rmin+width)*sqrt(2.0)
        else if (ioppbc .eq. 6 .or. ioppbc .eq. 7) then
c         Truncated octahedron
c         edge parameter is the distance of the truncating squares
c         from the center
c         Inscribed sphere is tangent to centers of opposite hexagons
          edge(1)=(rmin+width)*2.0/sqrt(3.0)
        else if (ioppbc .eq. 8) then
c         Hexagonal close packing
          edge(1)=2.0*(rmin+width)
        end if
      else if (ioppbc .gt. 0) then
c       Rectangular, cubic or hexagonal
        call extension(co,ih,nnh,1,n,cmin,cmax,c0,0,1,v)
        do k=1,3
          write (6,1004) xyz(k),c0(k),cmin(k),cmax(k)
          edge(k)=cmax(k)-cmin(k)+2.0*width
        end do
        if (ioppbc .eq. 2) then
          write (6,1000) (xyz(k),cmin(k),cmax(k),k=1,3)
        else if (ioppbc .eq. 1) then
          emax=0
          do k=1,3
            if (edge(k) .gt. emax) then
              emax=edge(k)
              kmax=k
            end if
          end do
          do k=1,3
            edge(k)=emax
          end do
          write (6,1001) xyz(kmax),cmin(kmax),cmax(kmax)
        else
          ix1=ixyzhex(1)
          write (6,1002) xyz(ix1),cmin(ix1),cmax(ix1),
     -     (xyz(ixyzhex(k)),k=2,3),
     -     (cmin(ixyzhex(k)),cmax(ixyzhex(k)),k=2,3)
           edge(2)=(amax1((cmax(ixyzhex(2))-cmin(ixyzhex(2))),
     -     (cmax(ixyzhex(3))-cmin(ixyzhex(3))))+2.0*width)*2.0/sqrt(3.0)
          edge(3)=edge(2)
        end if
        call shiftmol(co,n,c0,cn,-1.0)
      end if
      call prtcell(ioppbc,edge,edge_gen,0.0,volnew,nwnew,1)
      if (ireduce .eq. 1) then
        if (npbc .gt. 1) then
          call cellreduce(0,cn,ih,n,nnh,edge,edgen,edge_gen,ioppbc,npbc,
     -      cell,ncell,cellalt,ixyzhex,0.0,volnew,nwnew)
          call trnsfr(edge,edgen,3)
        end if
        call cellreduce(1,cn,ih,n,nnh,edge,edgen,edge_gen,ioppbc,npbc,
     -    cell,ncell,cellalt,ixyzhex,closlim,volnew,nwnew)
        call trnsfr(edge,edgen,3)
      end if
      return
1000  format(' Initial edges of the rectangle are based on the X, Y,',
     -  ' and Z extensions:',/,5x,3(1x,a1,f8.3,' - ',f8.3,5x))
1001  format(' Initial edge of the cube is based on the ',a1,
     -  ' extension:',f8.3,' - ',f8.3)
1002  format(' Initial length of the prism based the ',a,' extension:',
     -  ' (',f8.3,' - ',f8.3,')',/,' Initial edge of the hexagon ',
     -  'based on the ',a,' and ',a,' extensions:',
     -  /,' (',f8.3,' - ',f8.3,'), (',f8.3,' - ',f8.3,')')
1003  format(' Inscribed sphere radius=',f10.5)
1004  format(1x,a1,' c0=',f10.5,' cmin=',3f10.5,' cmax=',3f10.5)
      end
