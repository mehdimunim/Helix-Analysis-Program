      subroutine sortbondlist(ixres,listlen,indexbond,maxrsd,mxbonds)
      dimension ixres(maxrsd),indexbond(mxbonds)
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (6*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      parameter (IFILL6=MAX2D*MAX2D-6*MAXBONDS)
      common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),fill(IFILL6),
     -  i1(MAXBONDS),i2(MAXBONDS),i3(MAXBONDS),i4(MAXBONDS),
     -  newpair(2,MAXBONDS),ihbtores(MAXBONDS),nusepair(MAXBONDS),
     -  nhb_atot(MAXBONDS),nhb_rtot(MAXBONDS),a1(MAXBONDS),
     -  a2(MAXBONDS),fill_2(IFILL2)
      common /bondpairs/ ihbpair(2,MAXBONDS),ihb_pair_res(3,MAXBONDS)
      do i=1,listlen
        ii=indexbond(i)
        a1(i)=ixres(ihbpair(1,ii))*MAXBONDS+ixres(ihbpair(2,ii))
      end do
      call mrgsrt(6,indexbond,a1,listlen,i2,i3,i4,a2,listlen)
c      do i=1,listlen
c        ii=indexbond(i)
c        iii=ixres(ihbpair(1,ii))*MAXBONDS+ixres(ihbpair(2,ii))
c        write (6,5511) i,ii,(ixres(ihbpair(k,ii)),k=1,2),iii
c5511    format(i4,' ixb=',i3,' ir1,2=',2i4,' a=',i9)
c      end do
      return
      end
