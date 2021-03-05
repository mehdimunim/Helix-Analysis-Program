      subroutine averageres(nframe,res,ix,ires,maxframe,maxres,av,sd)
      dimension res(2,maxframe,maxres)
      real*8 sum,sum2
      sum=0.d0
      sum2=0.d0
      do iframe=1,nframe
        sum=sum+res(ix,iframe,ires)
        sum2=sum2+res(ix,iframe,ires)**2
      end do
      av=sum/nframe
      sd=dsqrt(dabs(sum2/nframe-(sum/nframe)**2))
      return
      end
