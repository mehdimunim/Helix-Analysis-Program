      subroutine residcorr(c1,c2,n,index,ncorr,nframe)
      dimension c1(3,n),c2(3,n),index(ncorr)
      parameter (MAXPHI=400,MAX2D=5000)
      parameter (IFILL4=MAXPHI*MAXPHI*MAXPHI-
     -  (2*MAX2D*MAX2D+17*MAX2D))
      real*8 trajcorr,cav1,cav2,cavs1,cavs2
      common /nnwork/ trajcorr(MAX2D,MAX2D),cav1(3,MAX2D),
     -  cav2(3,MAX2D),cavs1(MAX2D),cavs2(MAX2D),
     -  row(MAX2D),fill(IFILL4)
c     print *,'RESIDCORRR nframe,ncorr=',nframe,ncorr
      if (nframe .eq. 1) then
        call zeroitd(cav1,3*ncorr)
        call zeroitd(cav2,3*ncorr)
        call zeroitd(cavs1,ncorr)
        call zeroitd(cavs2,ncorr)
        do ir=1,ncorr
          call zeroitd(trajcorr(1,ir),ncorr)
        end do
      end if
      do ir=1,ncorr
        do k=1,3
          cav1(k,ir)=cav1(k,ir)+c1(k,index(ir))
          cav2(k,ir)=cav2(k,ir)+c2(k,index(ir))
          cavs1(ir)=cavs1(ir)+scprod(c1(k,index(ir)),c1(k,index(ir)))
          cavs2(ir)=cavs2(ir)+scprod(c2(k,index(ir)),c2(k,index(ir)))
        end do
        do jr=1,ir
          trajcorr(ir,jr)=trajcorr(ir,jr)+
     -      scprod(c1(1,index(ir)),c2(1,index(jr)))
          trajcorr(jr,ir)=trajcorr(ir,jr)
        end do
      end do
      return
      end
