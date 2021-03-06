      subroutine RMScalc(co,nats,dir,ip,fp,rms)
      implicit real*8(a-h,o-z)
      real*8 co(3,nats) !the input coordinates
c     Modified from the original to calculate the sd of the atom to axis dist
      integer nats              !number of atoms
      real*8 ip(3)              !the coordinates of the initial point
      real*8 fp(3)              !the coordinates of the final point
      real*8 dir(3)             !the axis vector

      real*8 tmp(3)

c     calculate rms deviations
      rms=0.0
      sum=0.0
      sum2=0.0
      do i=1,nats
         call dvdif (co(1,i),ip,tmp)
         call dvproj(dir,tmp,tmp)
         call dvsum(ip,tmp,fp)
         dd=ddistsq(co(1,i),fp)
         sum=sum+sqrt(dd)
         sum2=sum2+dd
      end do
      rms=sqrt(sum2/float(nats)-(sum/float(nats))**2)
      end
