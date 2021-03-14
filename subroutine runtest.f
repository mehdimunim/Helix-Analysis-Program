      subroutine runtest(nup,ndown,nrun,idecide)
c*****Computes test for correlation, type of correlation
      dimension nmncrt(20,20),nmxcrt(20,20)
c     Minimum critical values for correlation test:
      data nmncrt/20*0,
     -  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2,
     -  0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3,
     -  0, 0, 0, 0, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4,
     -  0, 0, 0, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5,
     -  0, 0, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6,
     -  0, 0, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6,
     -  0, 0, 2, 3, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7,
     -  0, 0, 2, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8,
     -  0, 0, 2, 3, 3, 4, 5, 5, 5, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9,
     -  0, 0, 2, 3, 4, 4, 5, 5, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 9,
     -  0, 2, 2, 3, 4, 4, 5, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9,10,10,
     -  0, 2, 2, 3, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 9,10,10,10,10,
     -  0, 2, 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 9, 9,10,10,10,11,11,
     -  0, 2, 3, 3, 4, 5, 6, 6, 7, 7, 8, 8, 9, 9,10,10,11,11,11,12,
     -  0, 2, 3, 4, 4, 5, 6, 6, 7, 8, 8, 9, 9,10,10,11,11,11,12,12,
     -  0, 2, 3, 4, 4, 5, 6, 7, 7, 8, 9, 9,10,10,11,11,11,12,12,13,
     -  0, 2, 3, 4, 5, 5, 6, 7, 8, 8, 9, 9,10,10,11,11,12,12,13,13,
     -  0, 2, 3, 4, 5, 6, 6, 7, 8, 8, 9,10,10,11,11,12,12,13,13,13,
     -  0, 2, 3, 4, 5, 6, 6, 7, 8, 9, 9,10,10,11,12,12,13,13,13,14/
c     Maximum critical values for correlation test:
      data nmxcrt/60*0,
     -  0, 0, 0, 0, 9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     -  0, 0, 0, 9,10,10,11,11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     -  0, 0, 0, 9,10,11,12,12,13,13,13,13, 0, 0, 0, 0, 0, 0, 0, 0,
     -  0, 0, 0, 0,11,12,12,12,14,14,14,14,15,15,15, 0, 0, 0, 0, 0,
     -  0, 0, 0, 0,11,12,13,13,14,15,15,16,16,16,16,17,17,17,17,17,
     -  0, 0, 0, 0, 0,13,14,14,15,16,16,16,17,17,18,18,18,18,18,18,
     -  0, 0, 0, 0, 0,13,14,15,16,16,17,17,18,18,18,19,19,19,20,20,
     -  0, 0, 0, 0, 0,13,14,15,16,17,17,18,19,19,19,20,20,20,21,21,
     -  0, 0, 0, 0, 0,13,14,16,16,17,18,19,19,20,20,21,21,21,22,22,
     -  0, 0, 0, 0, 0, 0,15,16,17,18,19,19,20,20,21,21,22,22,23,23,
     -  0, 0, 0, 0, 0, 0,15,16,17,18,19,20,20,21,22,23,23,23,23,24,
     -  0, 0, 0, 0, 0, 0,15,16,18,18,19,20,21,22,22,23,23,24,24,25,
     -  0, 0, 0, 0, 0, 0, 0,17,18,19,29,21,21,22,23,23,24,25,25,25,
     -  0, 0, 0, 0, 0, 0, 0,17,18,19,20,21,22,23,23,24,25,25,26,26,
     -  0, 0, 0, 0, 0, 0, 0,17,18,19,20,21,22,23,24,25,25,26,26,27,
     -  0, 0, 0, 0, 0, 0, 0,17,18,29,21,22,23,23,24,25,26,26,27,27,
     -  0, 0, 0, 0, 0, 0, 0,17,18,20,21,22,23,25,25,26,26,27,27,28/
      if (ndown .gt. 20 .or. nup .gt. 20) then
c       High
        idecide=5
      else if (ndown .eq. 0 .or. nup .eq. 0) then
c       Low
        idecide=1
      else if (nmxcrt(ndown,nup) .eq. 0 .or.
     -           nmncrt(ndown,nup) .eq. 0) then
c       Low
        idecide=1
      else if (nrun .lt. nmncrt(ndown,nup)) then
c       Correlated
        idecide=2
      else if (nrun .gt. nmxcrt(ndown,nup)) then
c       Correlated
        idecide=4
      else
c       Uncorrelated
        idecide=3
      end if
      return
      end
