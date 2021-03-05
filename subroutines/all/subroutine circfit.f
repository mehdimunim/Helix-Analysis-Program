      subroutine circfit(x,nats,dir,ip)
      implicit real*8 (a-h,o-z)
      real*8 x(3,nats)
c     Array of 3,2*MAXHX length
      real*8 xp(3,100)
      real*8 ax(6),sum(5),xp0(3)
      real*8 dir(3),ip(3),tmp(3)
      real*8 zero(3)
      integer nats
c     calculates the best-fit circle to a set of data and returns the
c     center in ip
c     x=the coordinates
c     dir(input) the direction vector
c     ip(output) the position vector
c     print *,'CIRCFIT nats,MAXHX=',nats,MAXHX

      data zero/0.,0.,0./
c     rotate the coordinates
c
c     read input
      call dvset(ax(1),ip)
      call dvset(ax(4),dir)

c      write(*,*) 'got ax:'
c      write(*,'(6(x,f8.3))') ax
c      write(*,*) 'got nats:',nats
c      write(*,*) 'got coords:'
c      do i=1,nats
c         write(*,'(3(x,f8.3))') (x(k,i),k=1,3)
c      end do

      call polar(dir(1),dir(2),dir(3),r,thet,fi)
c     write(*,*) 'Theta,phi',thet,fi
      do iat=1,nats
         call dvset(xp(1,iat),x(1,iat))
         call rotabout(xp(1,iat),zero,-thet,'z')
         call rotabout(xp(1,iat),zero,-fi,'y')
      end do
      call dvset(tmp,dir)
      do i=1,5
         sum(i)=0.
      end do

      do i=1,nats
         do j=i+1,nats
            sum(1)=sum(1)+(xp(1,i)-xp(1,j))**2
            sum(2)=sum(2)+(xp(2,i)-xp(2,j))**2
            sum(3)=sum(3)+(xp(1,i)-xp(1,j))*(xp(2,i)-xp(2,j))
            scrat=xp(1,i)**2+xp(2,i)**2-xp(1,j)**2-xp(2,j)**2
            sum(4)=sum(4)+scrat*(xp(1,i)-xp(1,j))
            sum(5)=sum(5)+scrat*(xp(2,i)-xp(2,j))
         end do
      end do
c     write(*,'(''sum'',5(x,f8.3))') sum

      xp0(1)=(sum(4)*sum(2)-sum(5)*sum(3))/
     -     (2.0*(sum(1)*sum(2)-sum(3)**2))
      xp0(2)=(sum(4)-2.0*xp0(1)*sum(1))/(2.0*sum(3))
      xp0(3)=xp(3,1)            !z coordinate of first atom
      call rotabout(xp0,zero,fi,'y')
      call rotabout(xp0,zero,thet,'z')
      call dvset(ax,xp0)

c     save output
      call dvset(ip,ax)
c     write(*,*) 'Output ip, dir'
c     write(*,'(6(x,f8.3))') ax
      end
