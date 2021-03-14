      subroutine roundlim(dmaxinp,div,ndiv)
c     For dmaxinp, find rounded div and ndiv such that div*ndiv ~ dmaxinp
      if (dmaxinp .gt. 1000000.0) then
        fac=10000.0
      else if (dmaxinp .gt. 100000.0) then
        fac=1000.0
      else if (dmaxinp .gt. 10000.0) then
        fac=100.0
      else if (dmaxinp .gt. 1000.0) then
        fac=10.0
      else
        fac=1.0
      end if
      dmax=dmaxinp/fac
      if (dmax .le. 0.05) then
        div=0.005
      else if (dmax .le. 0.1) then
        div=0.01
      else if (dmax .le. 0.2) then
        div=0.02
      else if (dmax .le. 0.5) then
        div=0.05
      else if (dmax .le. 1.0) then
        div=0.1
      else if (dmax .le. 2.0) then
        div=0.2
      else if (dmax .le. 5.0) then
        div=0.5
      else if (dmax .le. 8.0) then
        div=1.0
      else if (dmax .le. 10.0) then
        div=1.0
      else if (dmax .le. 12.0) then
        div=1.0
      else if (dmax .le. 15.0) then
        div=1.5
      else if (dmax .le. 20.0) then
        div=2.0
      else if (dmax .le. 25.0) then
        div=2.5
      else if (dmax .le. 50.0) then
        div=5.0
      else if (dmax .le. 100.0) then
        div=10.0
      else if (dmax .le. 200.0) then
        div=20.0
      else if (dmax .le. 250.0) then
        div=25.0
      else if (dmax .le. 400.0) then
        div=40.0
      else if (dmax .le. 500.0) then
        div=50.0
      else if (dmax .le. 600.0) then
        div=60.0
      else if (dmax .le. 750.0) then
        div=75.0
      else if (dmax .le. 800.0) then
        div=80.0
      else if (dmax .le. 1000.0) then
        div=100.0
      else
        div=200.0
      end if
      ndiv=(dmax+0.0001)/div
      if (dmax-ndiv*div .gt. dmax*0.01) ndiv=ndiv+1
      div=div*fac
      return
      end
