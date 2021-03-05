      subroutine plotbondcorr(ips,nbonds,yclab,lyclab,icontract,
     -  indexbond,title,temp,itemp,lablim,ncolcode,maxcolcode,
     -  ipspage,ipsclose)
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (6*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      dimension temp(MAX2D),itemp(MAX2D),indexbond(MAX2D)
      character*(*) yclab(nbonds)
      character*80 title
      common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
     -  ihbtores(MAXBONDS),nusepair(MAXBONDS),nhb_atot(MAXBONDS),
     -  nhb_rtot(MAXBONDS),ibond(MAXBONDS),a1(MAXBONDS),fill(IFILL2)
      dimension kc(1,1)
      real*8 dc(1,1)
      if (nbonds .lt. lablim) then
        nyclab=nbonds
      else
        nyclab=1
        lyclab=1
      end if
      inc=max0(1,500/nbonds)
      ixdel0=25
      iydel0=100
      ixdel=25
      iydel=200
      iytop=0
      do i=1,MAX2D
        temp(i)=i
        itemp(i)=i
      end do
      call plotmat(ips,kc,rmsd2d,dc,nbonds,nbonds,0,0,0,0,1,1,
     -  ixdel0,iydel0,0,iytop,0.0,1.0,ncolcode,maxcolcode,ixdel,iydel,
     -  inc,1.0,itemp,indexbond,indexbond,title,80,
     -  'Bond correlation distance measure matrix',40,1,' ',1,temp,
     -  yclab,nyclab,lyclab,1,MAX2D,1,MAX2D,MAX2D,ipspage,ipsclose)
      if (icontract .eq. 1) then
        iytop=iytop+15
        call pswrite(ips,ixdel0,iytop,'m',1)
        write (ips,*)
     -   '(Bond tracks contracted to residue-residue tracks) show'
      end if
      iydel=iydel-50
      ixcent=max0(0,(nbonds*inc-50*ncolcode)/2)
      call colcodeminmax(ips,25+ixcent,-iydel-5,0,ncolcode,
     -  maxcolcode,0.0,1.0)
      write (ips,1000) 'showpage'
      return
1000  format(a)
      end
