      subroutine opt2d(c,cnew,ih,n,nslt,nnh,edge,ioppbc,cell,ncell,
     -  ixyzhex,angle,rot,ftol,mintyp,i2dopt,iter)
      dimension c(3,n),cnew(3,n),ih(n),edge(3),cell(3,27),rot(3,3),
     -  angle(3),ixyzhex(3)
c     print *,'--- 2D optimization started'
      al=0.0
      angle(1)=al
      dl=touch(c,ih,nslt,nnh,edge,ioppbc,cell,ncell,ixyzhex,angle,rot,
     -  mintyp,i2dopt,cnew)
c     print *,'Initial touch=',dl,' i2dopt=',i2dopt
      ar=0.1
      angle(1)=ar
      dr=touch(c,ih,nslt,nnh,edge,ioppbc,cell,ncell,ixyzhex,angle,rot,
     -  mintyp,i2dopt,cnew)
      am=(al+ar)/2.0
      angle(1)=am
      dm=touch(c,ih,nslt,nnh,edge,ioppbc,cell,ncell,ixyzhex,angle,rot,
     -  mintyp,i2dopt,cnew)
      iter=0
c     print *,'al,am,ar=',al,am,ar
c     print *,'dl,dm,dr=',dl,dm,dr
      do while (abs(dl-dr) .gt. ftol)
        iter=iter+1
        if (dm .le. dl .and. dm .le. dr) then
c         Middle is the lowest
          if (dl .lt. dr) then
            ar=am
            am=(al+ar)/2.0
            dr=dm
          else
            al=am
            am=(al+ar)/2.0
            dl=dm
          end if
          angle(1)=am
          dm=touch(c,ih,nslt,nnh,edge,ioppbc,cell,ncell,ixyzhex,
     -      angle,rot,mintyp,i2dopt,cnew)
        else if (dl .le. dm .and. dl .le. dr) then
c         Left is the lowest
          d=-(am-al)
          ar=am
          dr=dm
          am=al
          dm=dl
          al=am+d
          angle(1)=al
          dl=touch(c,ih,nslt,nnh,edge,ioppbc,cell,ncell,ixyzhex,
     -      angle,rot,mintyp,i2dopt,cnew)
        else
c         Right is the lowest
          d=+(ar-am)
          al=am
          dl=dm
          am=ar
          dm=dr
          ar=am+d
          angle(1)=ar
          dr=touch(c,ih,nslt,nnh,edge,ioppbc,cell,ncell,ixyzhex,
     -      angle,rot,mintyp,i2dopt,cnew)
        end if
        if (mod(iter,05) .eq. 0) write (6,1000) iter,dm,abs(dl-dr),ftol
      end do
      return
1000  format(' 2D opt it',i5,' dm=',f9.5,' |dl-dr|=',f11.6,' tol=',f9.6)
      end
