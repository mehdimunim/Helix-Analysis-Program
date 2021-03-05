      subroutine dtqli(d,e,n,np,z,ierr)
      real*8 d(np),e(np),z(np,np)
      real*8 dd,r,g,s,c,p,f,b
      ierr=0
      if (n.gt.1) then
        do i=2,n
          e(i-1)=e(i)
        end do
        e(n)=0.d0
        do l=1,n
          iter=0
1         do m=l,n-1
            dd=abs(d(m))+abs(d(m+1))
c           if (abs(e(m))+dd.eq.dd) go to 2
c           MM 10/21/2004
            if (e(m) .eq. 0.0) go to 2
            if (dd .ne. 0.0) then
              if (abs(e(m))/dd .lt. 1.d-15) go to 2
            end if
          end do
          m=n
2         if (m. ne. l) then
c           if (iter .eq. 300) pause 'too many iterations'
            if (iter .eq. 300) then
              print *,'error: too many iterations in tqli'
              ierr=1
              return
            end if
            iter=iter+1
            g=(d(l+1)-d(l))/(2.d0*e(l))
            r=sqrt(g**2+1.d0)
            g=d(m)-d(l)+e(l)/(g+sign(r,g))
            s=1.d0
            c=1.d0
            p=0.d0
            do i=m-1,l,-1
              f=s*e(i)
              b=c*e(i)
              if(abs(f).ge.abs(g))then
                c=g/f
                r=dsqrt(c**2+1.d0)
                e(i+1)=f*r
                s=1.d0/r
                c=c*s
              else
                s=f/g
                r=dsqrt(s**2+1.d0)
                e(i+1)=g*r
                c=1.d0/r
                s=s*c
              endif
              g=d(i+1)-p
              r=(d(i)-g)*s+2.d0*c*b
              p=s*r
              d(i+1)=g+p
              g=c*r-b
              do k=1,n
                f=z(k,i+1)
                z(k,i+1)=s*z(k,i)+c*f
                z(k,i)=c*z(k,i)-s*f
              end do
            end do
            d(l)=d(l)-p
            e(l)=g
            e(m)=0.d0
            go to 1
          endif
        end do
      endif
      return
      end
