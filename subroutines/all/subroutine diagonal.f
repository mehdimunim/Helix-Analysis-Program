      subroutine diagonal(cn,ih,n,nnh,edge,ioppbc,cell,ncell,
     -  ixyzhex,closorgd,closest)
      dimension cn(3,n),ih(n),edge(3),cell(3,27),ixyzhex(3)
c     Orient molecule along the diagonal
      dimension rot(3,3)
      closorgd=distminimg(cn,ih,n,nnh,edge,ioppbc,cell,ncell,ixyzhex,
     -  inod,jnod)
c     Calculate rotation matrix
      sq2=sqrt(edge(1)**2+edge(2)**2)
      sq3=sqrt(edge(1)**2+edge(2)**2+edge(3)**2)
      rot(1,1)= edge(1)/sq3
      rot(1,2)=-edge(2)/sq3
      rot(1,3)=-edge(3)/sq3
      rot(2,1)= edge(2)/sq2
      rot(2,2)= edge(1)/sq2
      rot(2,3)=0.0
      rot(3,1)= edge(1)*edge(3)/(sq2*sq3)
      rot(3,2)=-edge(2)*edge(3)/(sq2*sq3)
      rot(3,3)= sq2/sq3
c     call chkort(rot)
      call rotate_c(cn,n,rot,cn,'DIAGONAL',8)
      closest=distminimg(cn,ih,n,nnh,edge,ioppbc,cell,ncell,ixyzhex,
     -  innd,jnnd)
      return
      end
