      subroutine kahn(co,nats,docircfit,dir,ip,fp,rms,
     -  iprint,message,MAXHX)
      implicit real*8(a-h,o-z)
c     ref. Computers in Chemistry Vol 13, No 3, pg 191, 1989
c     Approach: Construct a vector A from ca atom i to i-1. Construct B
c     from i to i+1.  Find V1, the vector which bisects A and B.  V1 is
c     perpendicuar to the helix axis (for a perfect mathematical helix).
c     Let i = i + 1 and repeat the procedure, finding V2.  As V1 and V2
c     are both perp. to the helical axis, their cross product gives the
c     helix direction.   Average over all possible tetrads of ca atoms,
c     or use a 3-D fitting method to give the direction based on the
c     points calculated to lie on the axis.
c
c     P1 and P2 are the vectors from the origin to CA 1 and CA 2.
c
c     The radius is a calculated as:
c     |dH|**2 - |p2-p1|**2
c     r = --------------------
c     2 * |(p1-p2) dot v2|
c     where d= (p2-p1) dot h
c
c     H1 and H2 are the position vectors
c     H1 = P1 + r*V1
c     H2 = P2 + r*V2
c

      real*8 co(3,MAXHX)
      integer nats,iprint
      real*8 ip(3)              !the coordinates of the initial point
      real*8 fp(3)              !the coordinates of the final point
      real*8 dir(3)             !the axis vector

      logical docircfit

      character*60 Message
      real*8 hsum(3)
c     Arrays of length MAXHX

      real*8 p1(3),p2(3),h(3),v1(3),v2(3),a(3),b(3),r,d
      real*8 tmp(3),p1mp2(3),p2mp1(3),dmag,ddot
c     Arrays of length 3,2*MAXHX
      real*8 h1(3,100)
      integer hcount

      if (iprint .gt. 3) write (6,7777) ((co(k,i),k=1,3),i=1,nats)
7777  format(' Axis routine input:',/,(3f10.4))

      do i=1,3
         hsum(i)=0.0
      end do
      hcount=1
      do i=2,nats-2
c     load p1 and p2
         call dvset(p1,co(1,i))
         call dvset(p2,co(1,i+1))
c     get vector v1
         call dvdif (p1,co(1,i-1),a)
         call dvdif (p1,co(1,i+1),b)
         call dvnorm(a)
         call dvnorm(b)
         call dvsum(a,b,v1)
         call dvnorm(v1)
c     get vector v2
         call dvdif (p2,co(1,i),a)
         call dvdif (p2,co(1,i+2),b)
         call dvnorm(a)
         call dvnorm(b)
         call dvsum(a,b,v2)
         call dvnorm(v2)

c     H=direction of axis
         call dcross(v1,v2,h)
         call dvnorm(h)
         call dvsum(h,hsum,hsum)

c     calculate radius
         call dvdif (p1,p2,p1mp2)
         call dvdif (p2,p1,p2mp1)
         d=ddot(p2mp1,h)        !rise/residue
         call dvmul(h,d,tmp)

c     kahn
         r=(dmag(tmp)**2 - dmag(p2mp1)**2)
     -        /(2.0*abs(ddot(p2mp1,v2)))

c     calculate the points on the axis
         call dvmul(V1,r,tmp)
         call dvsum(P1,tmp,H1(1,hcount))
         call dvmul(V2,r,tmp)
         call dvsum(P2,tmp,H1(1,hcount+1))
c     hcount=hcount+1
         hcount=hcount+2
      end do

      call dvset(dir,hsum)
      call dvnorm(dir)

      call parlsq(h1,hcount-1,docircfit,dir,ip,fp,rms,0)

      if (docircfit) then
         call circfit(co,nats,dir,ip)
      else
c     adjust ip to be next to the first alpha carbon, not the second
         call dvdif (co(1,1),ip,tmp)
         call dvproj(dir,tmp,tmp)
         call dvsum(ip,tmp,ip)
      end if

      call RMScalc(co,nats,dir,ip,fp,RMS)
      call writeout_h(dir,ip,fp,rms,message,iprint)
      return
      end
