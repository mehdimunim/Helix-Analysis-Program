      subroutine celldipole(c,n,nslt,index,nats,q,icharges,aw,iout,
     -  itraj)
      dimension c(3,n),index(n),q(n),aw(n)
c     Calculate the COM of the solute
c     Calculate the solute and the whole cell dipole moment when icharges >0
      parameter (MAXFRAMES=50000,MAXCOPY=600)
      common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
     -  increment,increment2,res(2,MAXFRAMES,MAXCOPY),
     -  x0(MAXCOPY),y0(MAXCOPY),nxselres,ixselres(MAXCOPY)
      real*8 dipole(3),com(3),awsum
      dimension dipslt(3),dipslv(3),diptot(3)
      call trajlimtest(nframe,MAXFRAMES)
      call zeroitd(dipole,3)
      call zeroitd(com,3)
      awsum=0.d0
      do ia=1,nats
        do k=1,3
          if (icharges .gt. 0)
     -      dipole(k)=dipole(k)+c(k,index(ia))*q(index(ia))
            com(k)=com(k)+c(k,index(ia))*aw(index(ia))
        end do
        awsum=awsum+aw(index(ia))
      end do
      do k=1,3
        if (icharges .gt. 0) dipslt(k)=dipole(k)
        com(k)=com(k)/awsum
      end do
      if (iout .gt. 0) write (iout,1003) com
      if (itraj .gt. 0) then
        res(2,nframe,14)=com(1)
        res(1,nframe,15)=com(2)
        res(2,nframe,15)=com(3)
      end if
      if (icharges .eq. 0) return
      dipsltabs=sqrt(dipslt(1)**2+dipslt(2)**2+dipslt(3)**2)
      call zeroitd(dipole,3)
      do ia=nslt+1,n
        do k=1,3
          dipole(k)=dipole(k)+c(k,ia)*q(ia)
        end do
      end do
      do k=1,3
        dipslv(k)=dipole(k)
      end do
      dipslvabs=sqrt(dipslv(1)**2+dipslv(2)**2+dipslv(3)**2)
      do ia=1,nslt
        do k=1,3
          dipole(k)=dipole(k)+c(k,index(ia))*q(index(ia))
        end do
      end do
      diptotabs=sqrt(diptot(1)**2+diptot(2)**2+diptot(3)**2)
      if (iout .ne. 0) then
        write (iout,1000) dipsltabs
        if (n .gt. nslt) write (iout,1001) dipslvabs,diptotabs
        write (iout,1002) 'solute ',dipslt
        if (n .gt. nslt) then
          write (iout,1002) 'solvent',dipslv
          write (iout,1002) 'cell   ',diptot
        end if
      end if
      if (itraj .gt. 0) then
        res(1,nframe,10)=diptotabs
        res(2,nframe,10)=dipsltabs
        res(1,nframe,11)=diptot(1)
        res(2,nframe,11)=diptot(2)
        res(1,nframe,12)=diptot(3)
        res(1,nframe,13)=dipslt(1)
        res(2,nframe,13)=dipslt(2)
        res(1,nframe,14)=dipslt(3)
      end if
      return
1000  format(' Dipole moment of the solute=',f10.3)
1001  format(' Dipole moment of the solvents=',f10.3,
     -  ' total=',f10.3,' au*A')
1002  format(' Dipole moment vector of ',a,'=',3f10.3,' au*A')
1003  format(' Center-of-mass=',3f10.4,' A')
      end
