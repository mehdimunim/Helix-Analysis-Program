      subroutine parlsq(co,n,docircfit,dir,ip,fp,rms,iprint)
c     calculates the axis of a helix from alpha carbon coordinates using
c     linear least squares regressions of atom number versus x,y, and z
c     coordinates.

      implicit real*8(a-h,o-z)!@Mehdi i-n are complex
      real*8 co(3,n)
      integer n,iprint          !n=number of atoms
      real*8 ip(3)              !the coordinates of the initial point
      real*8 fp(3)              !the coordinates of the final point
      real*8 dir(3)             !the axis vector
      logical docircfit

      character*60 message
c     print *,'PARLSQ MAXHX,n=',MAXHX,n
c      write (6,7777) ((co(k,i),k=1,3),i=1,n)
c7777  format(' PARLSQ input:',/,(3f10.4))
c      print *,'docircfit=',docircfit

      St=0.
      Srx=0.
      Sry=0.
      Srz=0.
      St_rx=0.
      St_ry=0.
      St_rz=0.
      St2=0.
      do i=1,n
         t= i-1
         St=St + t
         Srx=Srx + co(1,i)
         Sry=Sry + co(2,i)
         Srz=Srz + co(3,i)
         St_rx=St_rx + t * co(1,i)
         St_ry=St_ry + t * co(2,i)
         St_rz=St_rz + t * co(3,i)
         St2=St2 + t**2
      end do

      D=St2*n - St*St
      dir(1)=(n*St_rx - St*Srx)/D
      dir(2)=(n*St_ry - St*Sry)/D
      dir(3)=(n*St_rz - St*Srz)/D

      if (docircfit) then
c     use circular fit for ip
         call circfit(co,n,dir,ip)
      else
c     use least-squares initial point
         ip(1) =(St2*Srx - St*St_rx)/D
         ip(2) =(St2*Sry - St*St_ry)/D
         ip(3) =(St2*Srz - St*St_rz)/D
      end if

      call dvnorm(dir)
      call RMScalc(co,n,dir,ip,fp,RMS)
      message='parlsq'
      call writeout_h(dir,ip,fp,rms,message,iprint)
      end
