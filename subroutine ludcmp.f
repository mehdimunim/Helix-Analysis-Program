      subroutine ludcmp(a,n,np,indx,d)
      parameter (NMAX=100,TINY=1.0e-20)
      real a(np,np),vv(nmax)
      integer indx(n)
      data imax /0/
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) write(*,*) 'Oops, singular matrix!'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        if (j.gt.1) then
          do 14 i=1,j-1
            sum=a(i,j)
            if (i.gt.1) then
              do 13 k=1,i-1
                sum=sum-a(i,k)*a(k,j)
13            continue
              a(i,j)=sum
            end if
14        continue
        end if
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          if (j.gt.1) then
            do 15 k=1,j-1
              sum=sum-a(i,k)*a(k,j)
15          continue
            a(i,j)=sum
          end if
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          end if
16      continue
        if (j.ne.imax) then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        end if
        indx(j)=imax
        if (j.ne.n) then
          if (a(j,j).eq.0.)a(j,j)=tiny
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        end if
19    continue
      if (a(n,n).eq.0.) a(n,n)=tiny
      return
      end
