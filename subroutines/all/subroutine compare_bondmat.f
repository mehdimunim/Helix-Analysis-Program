      subroutine compare_bondmat(resnames,nrescol,itemp,bfac,nres,
     -  plotfile,lplotfile,mxrsd)
      dimension itemp(mxrsd),bfac(mxrsd)
      character*8 resnames(mxrsd)
      character(*) plotfile
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (6*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ rnb1(MAX2D,MAX2D),nng(MAX2D),rnb2(MAX2D,MAX2D),
     -  ihbtores(MAXBONDS),nusepair(MAXBONDS),nhb_atot(MAXBONDS),
     -  nhb_rtot(MAXBONDS),a1(MAXBONDS),a2(MAXBONDS),fill(IFILL2)
      common /colorinfo/ ncolcode,maxcolcode
      character*1 ans
      character*4 yclab(1)
      character*200 title1,title2,inpfile1,inpfile2
      dimension kc(1,1),xyval(1)
      real*8 dc(1,1),sum1,sum2
      data nyclab /1/,lyclab /1/
c     print *,'COMPARE_BONDMAT maxcolcode=',maxcolcode,' MAX2D=',MAX2D
      call read_bondmat(rnb1,41,'first',5,nres,inpfile1,linpfile1,
     -  itemp,MAX2D)
c     do iy=1,nres
c       write (77,8811) iy,(rnb1(ix,iy),ix=1,nres)
c     end do
      call read_bondmat(rnb2,41,'second',6,nres,inpfile2,linpfile2,
     -  itemp,MAX2D)
c     do iy=1,nres
c       write (78,8811) iy,(rnb2(ix,iy),ix=1,nres)
c     end do
c8811  format(' iy=',i5,/,(10f7.2))
      lplotfile=0
      write (6,1004)
      call openfile(42,0,'residue bond matrix difference',30,'new',
     -  plotfile,lplotfile,notfnd,0,1,1,0,0)
      call quiz(ans,idifftyp,' ',' ',0,'bond-difference plot type',25,
     -  0,5,6,000)
      rnbmin=99999.9
      rnbmax=-rnbmin
      sum1=0.d0
      sum2=0.d0
      isum1=0
      isum2=0
      call zeroit(bfac,nres)
      do ir=1,nres
        do in=1,nres
          sum1=sum1+rnb1(ir,in)
          sum2=sum2+rnb2(ir,in)
          if (rnb1(ir,in) .gt. 0.0 .and. rnb2(ir,in) .eq. 0.0)
     -      isum1=isum1+1
          if (rnb1(ir,in) .eq. 0.0 .and. rnb2(ir,in) .gt. 0.0)
     -      isum2=isum1+1
          diff=rnb1(ir,in)-rnb2(ir,in)
          write (42,1000) ir,resnames(ir)(1:nrescol),
     -      in,resnames(in)(1:nrescol),diff
          if (diff .lt. rnbmin) rnbmin=diff
          if (diff .gt. rnbmax) rnbmax=diff
          if (idifftyp .eq. 1) then
            rnb1(ir,in)=diff
          else if (idifftyp .eq. 2) then
            rnb1(ir,in)=abs(diff)
          else if (idifftyp .eq. 3) then
            if (rnb1(ir,in) .eq. 0.0 .and. rnb2(ir,in) .gt. 0.0) then
              rnb1(ir,in)=+1.0
            else if (rnb1(ir,in) .gt. 0.0 .and.rnb2(ir,in) .eq. 0.0)then
              rnb1(ir,in)=-1.0
            else
              rnb1(ir,in)=0.0
            end if
          end if
          bfac(ir)=bfac(ir)+rnb1(ir,in)
        end do
c       bfac(ir)=bfac(ir)/float(nres)
      end do
      if (idifftyp .eq. 2) then
        rnbmax=amax1(abs(rnbmin),abs(rnbmax))
        rnbmin=0.0
      else if (idifftyp .eq. 3) then
        rnbmax=1.0
        rnbmin=-1.0
      end if
      write (6,1003) inpfile1(1:linpfile1),sum1,isum1
      write (6,1003) inpfile2(1:linpfile2),sum2,isum2
      write (6,1002) rnbmin,rnbmax
      if (idifftyp .ne. 3) then
        call getreal('Minimum difference to show',26,rnbmin,rnbmin,1,0)
        call getreal('Maximum difference to show',26,rnbmax,rnbmax,1,0)
      end if
      close (42)
      ips=43
      plotfile(lplotfile+1:lplotfile+3)='.ps'
      lplotfile=lplotfile+3
      call openfile(ips,0,'residue bond matrix difference',30,'new',
     -  plotfile,lplotfile,notfnd,0,1,1,0,0)
      write (6,1001) plotfile(1:lplotfile)
      title2='Average distance matrix difference plot'
      ltitle2=39
      if (idifftyp .eq. 2) then
        title2(ltitle2+1:ltitle2+18)=' - absolute values'
        ltitle2=ltitle2+18
      else if (idifftyp .eq. 3) then
        title2(ltitle2+1:ltitle2+27)=' - signs of the differences'
        ltitle2=ltitle2+27
      end if
      ltitle1=linpfile1+linpfile2+3
      title1(1:ltitle1)=
     -  inpfile1(1:linpfile1)//' - '//inpfile2(1:linpfile2)
      call openps(ips,500.0,500.0,' ',1,' ',1,plotfile,0,
     -  plotfile,0,1,ipspage)
      call getint('Number of colors to use',23,5,1,maxcolcode,ncolcode,
     -  0)
      scalefac=amin1(1.0,500.0/float(nres))
      ixdel=25
      iydel=115
      iytop=0
      incinp=max0(1,500/nres)
      call indexit(itemp,1,mxrsd,0)
      call plotmat(ips,kc,rnb1,dc,nres,nres,0,0,0,0,1,0,0,iydel,00,
     -  iytop,rnbmin,rnbmax,ncolcode,maxcolcode,ixdel,iydel,incinp,
     -  scalefac,itemp,itemp,itemp,title1,ltitle1,title2,ltitle2,0,' ',
     -  0,xyval,yclab,nyclab,lyclab,1,MAX2D,1,mxrsd,mxrsd,ipspage,0)
      ixd=ixdel
      if (ncolcode .le. 5) inxd=ixd+60
      call colcodeminmax(ips,ixd,-60,0,ncolcode,maxcolcode,rnbmin,
     -  rnbmax)
      return
1000  format(' Res #',i6,'(',a,') - Res #',i6,'(',a,'): <nb1>-<nb2>=',
     -  f6.3)
1001  format(' Residue bond difference matrix plots are on file ',a)
1002  format(' Range of differences: [',f6.3,',',f6.3,'] A')
1003  format(' Sum of bonds in ',a,'=',f10.3,' sum of (+) changes=',i6)
1004  format(' Comparing two residue-residue bond frequency matrices',/,
     -  ' written by the bond-tracking option of Simulaid')
      end
