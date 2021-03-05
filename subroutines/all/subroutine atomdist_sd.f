      subroutine atomdist_sd(c,n,ldist,ndist,nframe)
      dimension c(3,n),ldist(ndist)
      parameter (MAXPHI=400,MAX2D=5000)
      parameter (IFILL9=MAXPHI*MAXPHI*MAXPHI-
     -  (2*MAX2D*MAX2D+11*MAX2D))
      real*8 trajdist,cav,cavs
      common /nnwork/ trajdist(MAX2D,MAX2D),cav(3,MAX2D),
     -  cavsng(3,MAX2D),cavs(MAX2D),fill(IFILL9)
      real*8 d2
c     print *,'ATOMDIST_SD nframe,ndist=',nframe,ndist
      if (nframe .eq. 1) then
        call zeroitd(cav,3*ndist)
        call zeroitd(cavs,ndist)
        do ir=1,ndist
          call zeroitd(trajdist(1,ir),ndist)
        end do
      end if
      do ir=1,ndist
        do k=1,3
          cav(k,ir)=cav(k,ir)+c(k,ldist(ir))
          cavs(ir)=cavs(ir)+scprod(c(1,ldist(ir)),c(1,ldist(ir)))
        end do
        do jr=1,ir
          d2=dist2(c(1,ldist(ir)),c(1,ldist(jr)))
          trajdist(ir,jr)=trajdist(ir,jr)+d2
          trajdist(jr,ir)=trajdist(jr,ir)+dsqrt(d2)
        end do
      end do
      return
      end
