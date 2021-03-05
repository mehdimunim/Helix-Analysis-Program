      subroutine print_rrdist(itypavg,nframe,irefres1,irefres2,inegres1,
     -  inegres2,listrefres,nrefres,listnegres,nnegres,nrescol,iwrrdr,
     -  iwrrdc,ips,ipspage,resnames,inpfile,namleni,itemp,mxrsd)
      dimension itemp(mxrsd),listrefres(mxrsd),listnegres(mxrsd)
      character*8 resnames(mxrsd)
      character*(*) inpfile
      common /colorinfo/ ncolcode,maxcolcode
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (6*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ sdm(MAX2D,MAX2D),nng(MAX2D),rmsd(MAX2D,MAX2D),
     -  ihbtores(MAXBONDS),nusepair(MAXBONDS),nhb_atot(MAXBONDS),
     -  nhb_rtot(MAXBONDS),a1(MAXBONDS),a2(MAXBONDS),fill(IFILL2)
      dimension kc(1,1),xyval(1)
      character*4 yclab(1)
      character*14 distlab
      character*80 title2
      real*8 dc(1,1)
      data nyclab /1/,lyclab /1/
      if (ips .eq. 0) return
      if (itypavg .eq. 1) then
        iw=iwrrdr
        distlab='representative'
        ldistlab=14
      else
        iw=iwrrdc
        distlab='closest'
        ldistlab=7
      end if
      title2='Maximum '//distlab(1:ldistlab)//
     -  ' atom based distance to print'
      ltitle2=ldistlab+37
      call getreal(title2,ltitle2,999999.0,rprtmax,1,0)
      title2='Average residue-residue distances based on '//
     -  distlab(1:ldistlab)//' atom distances'
      ltitle2=ldistlab+62
      write (iw,1001) title2(1:ltitle2)
      do irrr=1,nrefres
        irr=listrefres(irrr)
        do inrr=1,nnegres
          inr=listnegres(inrr)
          rmsd(irr-irefres1+1,inr-inegres1+1)=
     -      rmsd(irr-irefres1+1,inr-inegres1+1)/float(nframe)
          sdm(irr-irefres1+1,inr-inegres1+1)=
     -    sqrt(sdm(irr-irefres1+1,inr-inegres1+1)/float(nframe)-
     -    rmsd(irr-irefres1+1,inr-inegres1+1)**2)
          if (rmsd(irr-irefres1+1,inr-inegres1+1) .le. rprtmax)
     -      write (iw,1000) irr,resnames(irr)(1:nrescol),inr,
     -        resnames(inr)(1:nrescol),
     -        rmsd(irr-irefres1+1,inr-inegres1+1),
     -        sdm(irr-irefres1+1,inr-inegres1+1)
        end do
      end do
      call getreal(
     -  'Upper limit of the distance range to color code in the matrix',
     -   61,999999.0,rcmax,1,0)
      rcmin=0.0
      call getint('Number of colors to use',23,5,1,maxcolcode,ncolcode,
     -  0)
      nresx=irefres2-irefres1+1
      nresy=inegres2-inegres1+1
c      do in=1,nresy
c        write (iw,9855) in+inegres1-1,(rmsd(ir,in),ir=1,nresx)
c9855    format(' in=',i4,200f6.2)
c      end do
      scalefac=amin1(1.0,500.0/float(max0(nresx,nresy)))
      ixdel=25
      iydel=115
      iytop=0
      incinp=max0(1,500/max0(nresx,nresy))
      call indexit(itemp,1,mxrsd,0)
      call plotmat(ips,kc,rmsd,dc,nresx,nresy,irefres1-1,inegres1-1,
     -  irefres1-1,inegres1-1,1,0,0,iydel,00,iytop,rcmin,rcmax,ncolcode,
     -  maxcolcode,ixdel,iydel,incinp,scalefac,itemp,itemp,itemp,
     -  inpfile,namleni,title2,ltitle2,0,' ',0,xyval,yclab,nyclab,
     -  lyclab,1,MAX2D,1,mxrsd,mxrsd,ipspage,0)
      ixd=ixdel
      if (ncolcode .le. 5) ixd=ixd+60
c     print *,'RCMIN,RCMAX=',rcmin,rcmax
      call colcodeminmax(ips,ixd,-60,0,ncolcode,maxcolcode,rcmin,rcmax)
      return
1000  format(' Res #',i6,'(',a,') - Res #',i6,'(',a,'): <d>=',f8.4,
     -  ' A sd=',f8.4)
1001  format(/,1x,a,' averages over the trajectory',/)
      end
