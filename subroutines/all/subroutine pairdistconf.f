      subroutine pairdistconf(c,n,rmsdmin,rmsdmax)
      dimension c(3,n)
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (6*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
     -  ihbtores(MAXBONDS),nusepair(MAXBONDS),nhb_atot(MAXBONDS),
     -  nhb_rtot(MAXBONDS),a1(MAXBONDS),a2(MAXBONDS),fill(IFILL2)
      rmsd2d(1,1)=0.0
      rmsdmax=0.0
      rmsdmin=99999.0
      do ia=2,n
        do ja=1,ia-1
          d2=sqrt(dist2(c(1,ia),c(1,ja)))
          rmsd2d(ia,ja)=d2
          rmsd2d(ja,ia)=d2
          if (d2 .lt. rmsdmin) rmsdmin=d2
          if (d2 .gt. rmsdmax) rmsdmax=d2
        end do
        rmsd2d(ia,ia)=0.0
      end do
      return
      end
