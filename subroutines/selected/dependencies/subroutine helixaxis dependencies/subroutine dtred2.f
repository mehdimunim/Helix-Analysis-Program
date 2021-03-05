      subroutine dtred2(a,n,np,d,e)
      real*8 a(np,np),d(np),e(np)
      real*8 scale,h,f,g
      if (n .gt. 1) then
        do i=n,2,-1
          l=i-1
          h=0.d0
          scale=0.d0
          if (l. gt. 1) then
            do k=1,l
              scale=scale+abs(a(i,k))
            end do
            if(scale.eq.0.)then
              e(i)=a(i,l)
            else
              do k=1,l
                a(i,k)=a(i,k)/scale
                h=h+a(i,k)**2
              end do
              f=a(i,l)
              g=-sign(sqrt(h),f)
              e(i)=scale*g
              h=h-f*g
              a(i,l)=f-g
              f=0.
              do j=1,l
                a(j,i)=a(i,j)/h
                g=0.
                do k=1,j
                  g=g+a(j,k)*a(i,k)
                end do
                if (l .gt. j) then
                  do k=j+1,l
                    g=g+a(k,j)*a(i,k)
                  end do
                end if
                e(j)=g/h
                f=f+e(j)*a(i,j)
              end do
              hh=f/(h+h)
              do j=1,l
                f=a(i,j)
                g=e(j)-hh*f
                e(j)=g
                do k=1,j
                  a(j,k)=a(j,k)-f*e(k)-g*a(i,k)
                end do
              end do
            endif
          else
            e(i)=a(i,l)
          end if
          d(i)=h
        end do
      endif
      d(1)=0.d0
      e(1)=0.d0
      do i=1,n
        l=i-1
        if (d(i) .ne. 0.d0) then
          do j=1,l
            g=0.
            do k=1,l
              g=g+a(i,k)*a(k,j)
            end do
            do k=1,l
              a(k,j)=a(k,j)-g*a(k,i)
            end do
          end do
        endif
        d(i)=a(i,i)
        a(i,i)=1.d0
        if (l .ge. 1) then
          do j=1,l
            a(i,j)=0.d0
            a(j,i)=0.d0
          end do
        endif
      end do
      return
      end
