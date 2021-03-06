      subroutine ramachandran_init(n,ixres)
      parameter (MAXREC=200000,MAXRSD=70000,MAXNEIG=70,MAXPHI=400)
      parameter (IFILL1=MAXPHI*MAXPHI*MAXPHI-(7+2*MAXNEIG)*MAXREC)
      parameter (IFILL5=(MAXNEIG+6)*MAXREC+IFILL1-44*MAXRSD)
      common /nnwork/ ineig(MAXNEIG,MAXREC),nneig(MAXREC),
     -  nprossacc(6,6,MAXRSD),issprossacc(5,MAXRSD),
     -  ixypross(2,MAXRSD),isspross(MAXRSD),fill(IFILL5)
      dimension ixres(n)
      nres=ixres(n)
c     print *,'RAMA init n,nres=',n,nres
      call zeroiti(nprossacc,0,36*nres)
      call zeroiti(issprossacc,0,5*nres)
      return
      end
