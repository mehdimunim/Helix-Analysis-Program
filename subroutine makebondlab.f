      subroutine makebondlab(ilab1,ilab2,inc,ibondlab,bondlab,lbondlab,
     -  irrix,ixresno,ixres,index,resnames,atnames,maxbondlab,mxbonds,
     -  maxrsd,maxrec)
      dimension ixresno(maxrsd),irrix(maxrsd),ixres(maxrec),
     -  index(mxbonds)
      character*8 atnames(maxrec),resnames(maxrsd)
      character*37 bondlab(maxbondlab)
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (6*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
     -  ihbtores(MAXBONDS),nusepair(MAXBONDS),nhb_atot(MAXBONDS),
     -  nhb_rtot(MAXBONDS),ibond(MAXBONDS),a1(MAXBONDS),fill(IFILL2)
      common /bondpairs/ ihbpair(2,MAXBONDS),ihb_pair_res(3,MAXBONDS)
c     print *,'MAKEBONDLAB ilab1,ilab2,maxbondlab,ilab=',
c    -  ilab1,ilab2,maxbondlab,ilab
      lbondlab=0
      do ilab=ilab1,ilab2
        call blankout(bondlab(ilab-inc),1,37)
        if (ibondlab .eq. 2) then
c         res-aggregated plot
          irr=irrix(ilab)
          ixr1=ihb_pair_res(1,irr)
          ixr2=ihb_pair_res(2,irr)
          ir1=ixresno(ixr1)
          ir2=ixresno(ixr2)
c         write (6,9871) ilab,irr,ixr1,ixr2,ir1,ir2
c9871     format(' ILAB,IRR=',2i4,' IXR1,IXR2=',2i4,' IR1,IR2=',2i4)
          call lastchar(resnames(ixr1),lc1,8)
          call lastchar(resnames(ixr2),lc2,8)
          lbondlabi=11+lc1+lc2
          write (bondlab(ilab-inc)(1:lbondlabi),1002)
     -      resnames(ixr1)(1:lc1),ir1,resnames(ixr2)(1:lc2),ir2
        else
c         Single bond plot
          ia1=ihbpair(1,index(ilab))
          ia2=ihbpair(2,index(ilab))
          ixr1=ixres(ia1)
          ixr2=ixres(ia2)
          ir1=ixresno(ixr1)
          ir2=ixresno(ixr2)
          call lastchar(resnames(ixr1),lc1,8)
          call lastchar(resnames(ixr2),lc2,8)
          lbondlabi=21+lc1+lc2
          write (bondlab(ilab-inc)(1:lbondlabi),1001) atnames(ia1)(1:4),
     -      resnames(ixr1)(1:lc1),ir1,atnames(ia2)(1:4),
     -      resnames(ixr2)(1:lc2),ir2
        end if
        if (lbondlab .lt. lbondlabi) lbondlab=lbondlabi
      end do
      return
1001  format(a,1x,a,i4,' - ',a,1x,a,i4)
1002  format(a,i4,' - ',a,i4)
      end
