      subroutine roundlimint(nmaxinp,idiv,ndiv)
c     For nmaxinp, find rounded idiv and ndiv such that idiv*ndiv ~ nmaxinp
      if (nmaxinp .gt. 1000000) then
        ifac=10000
      else if (nmaxinp .gt. 100000) then
        ifac=1000
      else if (nmaxinp .gt. 10000) then
        ifac=100
      else if (nmaxinp .gt. 1000) then
        ifac=10
      else
        ifac=1
      end if
      nmax=nmaxinp/ifac
      if (nmax .le. 10) then
        idiv=1
      else if (nmax .le. 50) then
        idiv=5
      else if (nmax .le. 100) then
        idiv=10
      else if (nmax .le. 200) then
        idiv=20
      else if (nmax .le. 250) then
        idiv=25
      else if (nmax .le. 300) then
        idiv=30
      else if (nmax .le. 400) then
        idiv=40
      else if (nmax .le. 500) then
        idiv=50
      else if (nmax .le. 600) then
        idiv=60
      else if (nmax .le. 750) then
        idiv=75
      else if (nmax .le. 800) then
        idiv=80
      else if (nmax .le. 1000) then
        idiv=100
      else
        idiv=200
      end if
c     print *,'ROUND idiv,ifac,nmaxinp=',idiv,ifac,nmaxinp
      idiv=idiv*ifac
      ndiv=(nmaxinp-1)/idiv+1
c     print *,'ROUND2 idiv,ndiv=',idiv,ndiv
      return
      end
