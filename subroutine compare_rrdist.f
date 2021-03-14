      subroutine compare_rrdist(resnames,nrescol,itemp,bfac,irefres1,
     -  irefres2,plotfile,lplotfile,mxrsd)
      dimension itemp(mxrsd),bfac(mxrsd)
      character*8 resnames(mxrsd)
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (6*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ rr1(MAX2D,MAX2D),nng(MAX2D),rr2(MAX2D,MAX2D),
     -  ihbtores(MAXBONDS),nusepair(MAXBONDS),nhb_atot(MAXBONDS),
     -  nhb_rtot(MAXBONDS),a1(MAXBONDS),a2(MAXBONDS),fill(IFILL2)
      common /colorinfo/ ncolcode,maxcolcode
      character*4 yclab(1)
      character*200 plotfile,title2,inpfile1,inpfile2
      character*80 title1
      dimension kc(1,1),xyval(1)
      real*8 dc(1,1)
      data nyclab /1/,lyclab /1/
c     print *,'COMPARE_RRDIST maxcolcode,mxrsd=',maxcolcode,mxrsd
      ifail=1
      do while (ifail .eq. 1)
        call read_rrdist(rr1,41,'first',5,irefres1,irefres2,inegres1,
     -    inegres2,inpfile1,linpfile1,ifail1,MAX2D)
        call read_rrdist(rr2,41,'second',6,irefres11,irefres22,
     -    inegres11,inegres22,inpfile2,linpfile2,ifail2,MAX2D)
        if (irefres1 .ne. irefres11 .or. irefres2 .ne. irefres22 .or.
     -      inegres1 .ne. inegres11 .or. inegres2 .ne. inegres22) then
          print *,'Residue ranges differ - can not run the comparison'
        else
          ifail=ifail1+ifail2
        end if
        if (ifail .gt. 0) write (6,1005)
      end do
      lplotfile=0
      call openfile(42,0,'average distance matrix difference',34,'new',
     -  plotfile,lplotfile,notfnd,0,1,1,0,0)
      rdmin=99999.0
      rdmax=-rdmin
      call askyn('Do you want absolute values in the difference plot',
     -  50,1,-1,iabsval,0,0)
      call getreal('Maximum distance to show',24,10000.0,dmaxshow,1,132)
      call zeroit(bfac,mxrsd)
      ncount_tot=0
      do irr=irefres1,irefres2
        ir=irr-irefres1+1
        ncount=0
        do inr=inegres1,inegres2
          in=inr-inegres1+1
          diffrn=rr1(ir,in)-rr2(ir,in)
          write (42,1000) irr,resnames(irr)(1:nrescol),
     -      inr,resnames(inr)(1:nrescol),rr1(ir,in),rr2(ir,in),diffrn
          if (iabsval .eq. 1) diffrn=abs(diffrn)
          if (diffrn .lt. rdmin) rdmin=diffrn
          if (diffrn .gt. rdmax) rdmax=diffrn
          if (rr1(ir,in) .gt. dmaxshow .and.
     -        rr2(ir,in) .gt. dmaxshow) then
            rr1(ir,in)=0.0
          else
            rr1(ir,in)=diffrn
            bfac(irr)=bfac(irr)+abs(rr1(ir,in))
            ncount=ncount+1
          end if
        end do
        ncount_tot=ncount_tot+ncount
        write (42,1002) ' N',dmaxshow,ncount
        if (ncount .eq. 0) ncount=1
        bfac(irr)=bfac(irr)/float(ncount)
      end do
      write (42,1002) ' Total n',dmaxshow,ncount_tot
      write (42,1003) dmaxshow
      do irr=irefres1,irefres2
        ir=irr-irefres1+1
        do inr=inegres1,inegres2
          in=inr-inegres1+1
          if (rr1(ir,in) .ne. 0.0) write (42,1004) irr,
     -      resnames(irr)(1:nrescol),inr,resnames(inr)(1:nrescol),
     -      rr1(ir,in),rr2(ir,in)
        end do
      end do
      if (iabsval .eq. 1) rdmin=0.0
      close (42)
      ips=43
      plotfile(lplotfile+1:lplotfile+3)='.ps'
      lplotfile=lplotfile+3
      call openfile(ips,0,'average distance matrix difference',34,'new',
     -  plotfile,lplotfile,notfnd,0,1,1,0,0)
      write (6,1001) plotfile(1:lplotfile),rdmin,rdmax
      title2='Average distance matrix difference plot'
      ltitle2=39
      ltitle1=linpfile1+linpfile2+3
      title1(1:ltitle1)=
     -  inpfile1(1:linpfile1)//' - '//inpfile2(1:linpfile2)
      call openps(ips,500.0,500.0,' ',1,' ',1,plotfile,0,
     -  plotfile,0,1,ipspage)
      call getint('Number of colors to use',23,5,1,maxcolcode,ncolcode,
     -  0)
      nresx=irefres2-irefres1+1
      nresy=inegres2-inegres1+1
      scalefac=amin1(1.0,500.0/float(max0(nresx,nresy)))
      ixdel=25
      iydel=115
      iytop=0
      incinp=max0(1,500/max0(nresx,nresy))
      call indexit(itemp,1,mxrsd,0)
      call plotmat(ips,kc,rr1,dc,nresx,nresy,irefres1-1,inegres1-1,0,0,
     -  1,0,0,iydel,00,iytop,rdmin,rdmax,ncolcode,maxcolcode,ixdel,
     -  iydel,incinp,scalefac,itemp,itemp,itemp,title1,ltitle1,title2,
     -  ltitle2,0,' ',0,xyval,yclab,nyclab,lyclab,1,MAX2D,1,mxrsd,
     -  mxrsd,ipspage,0)
      ixd=ixdel
      if (ncolcode .le. 5) inxd=ixd+60
      call colcodeminmax(ips,ixd,-60,0,ncolcode,maxcolcode,rdmin,rdmax)
      return
1000  format(' Res #',i6,'(',a,') - Res #',i6,'(',a,'):<d1>=',f6.2,
     -  ' <d2>=',f6.2,' <d1>-<d2>=',f7.3)
1001  format(' Residue distance difference matrix plots are on file ',a,
     -  /,' Range of differences: [',f10.5,',',f10.5,'] A')
1002  format(a,'umber of pairs within ',f5.1,' A:',i5)
1003  format(' Distance changes for pairs closer than ',f5.2,' A')
1004  format(' Res #',i6,'(',a,') - Res #',i6,'(',a,'):<d1>-<d2>=',f6.2,
     -  ' <d2>=',f6.2)
1005  format(' NOTE: input file should be the .rsd file written by the',
     -  /,'R<E>sidue distance list option of distance analysis')
      end
