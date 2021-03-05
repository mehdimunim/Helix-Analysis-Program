      subroutine res_res_bond(nres,nbfound,indexbond,nhbdist,iresno,
     -  ixres,ifres,isegno,irrix,it2,it3,title,ltitle,bondname,
     -  lbondname,line,index,irescol1,irescol2,iresshift,iout,iplot,
     -  ipspage,ncolcode,maxcolcode,mxbonds,maxrsd,maxat)
      dimension nhbdist(mxbonds),iresno(maxat),ixres(maxat),
     -  indexbond(mxbonds),ifres(maxat),isegno(maxat),index(maxat),
     -  irrix(maxat),it2(maxat),it3(maxat)
      character*(*) title,bondname
      character* 132 line(maxat)
      character*200 trajnam,trajnam2
      common /trajname/ trajnam,trajnam2,ltrajnam,ltrajnam2,ifirsttraj,
     -  ifirsttraj2,ilasttraj,ilasttraj2,incrementtraj,incrementtraj2
      parameter (MAXFRAMES=50000,MAXCOPY=600,MAXITEMS=2*MAXCOPY-2)
      common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
     -  increment,increment2,ires(MAXITEMS,MAXFRAMES),
     -  scres(2,MAXFRAMES),x0(MAXCOPY),y0(MAXCOPY),nxselres,
     -  ixselres(MAXCOPY)
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (6*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
     -  ihbtores(MAXBONDS),nusepair(MAXBONDS),nhb_atot(MAXBONDS),
     -  nhb_rtot(MAXBONDS),ibond(MAXBONDS),a1(MAXBONDS),fill(IFILL2)
      common /bondpairs/ ihbpair(2,MAXBONDS),ihb_pair_res(3,MAXBONDS)
      dimension ixshuffle(MAX2D)
      real*8 dc(1,1)
      character*1 ans
      character*13 yclab(40)
      character*80 question
c     print *,'RES_RES nframe,NBFOUND,MAXBONDS=',nframe,nbfound,MAXBONDS
      call indexit(ixshuffle,1,MAX2D,0)
      call zeroiti(irrix,0,nres)
      call zeroiti(it2,0,nres)
      call zeroiti(it3,0,nres)
      do i=1,nbfound
        ihb=indexbond(i)
c       write (iout,9876) i,ihb,(ihbpair(k,ihb),k=1,2),
c    -    (ixres(ihbpair(k,ihb)),k=1,2)
c9876   format(' I,IHB=',2i5,' IHBPAIR=',2i6,' IXRES=',2i6)
        irrix(ixres(ihbpair(1,ihb)))=1
        irrix(ixres(ihbpair(2,ihb)))=1
      end do
c     irrix(ir) is one if residue # ir is in a bond listed
c     it2 points from the original list to the short list
c     it3 points from the short residue list to the original list
      nreshb=0
      do ir=1,nres
        if (irrix(ir) .gt. 0) then
          nreshb=nreshb+1
          it2(ir)=nreshb
          it3(nreshb)=ir
        end if
      end do
      write (iout,1021)
      write (6,1018) bondname(1:lbondname),nreshb
      write (question,2000) bondname(1:lbondname)
      lq=45+lbondname
      call askyn(question,lq,1,+1,ireduceplot,73,0)
      call quiz(ans,iresres,'c',' ',0,'treatment of contacts',21,
     -  0,5,6,70)
      write (iout,1018) bondname(1:lbondname),nreshb
      write (iout,1000) bondname(1:lbondname),nbfound,nframe
      if (ireduceplot .eq. 0) then
        call indexit(it2,1,nres,0)
        call indexit(it3,1,nres,0)
        nreshb=nres
      end if
      if (ireduceplot .gt. 0 .or. iresshift .gt. 0) then
        write (iout,1001) bondname(1:lbondname),(ir,
     -    iresno(ifres(it3(ir))),isegno(ifres(it3(ir))),
     -    ixres(ifres(it3(ir))),ir=1,nreshb)
        if (iresshift .gt. 0) write (iout,1019)
        write (iout,1014)(it3(ir),ir=1,nreshb)
      end if
      nhbresmax=0
      do ir=1,nreshb
        do jr=1,nreshb
          ing(ir,jr)=0
          rmsd2d(ir,jr)=0.0
        end do
      end do
c     Calculate cumulative res-res bond percentage
      do i=1,nbfound
        ihb=indexbond(i)
        ir1=it2(ixres(ihbpair(1,ihb)))
        ir2=it2(ixres(ihbpair(2,ihb)))
        ing(ir1,ir2)=ing(ir1,ir2)+1
        if (nhbresmax .lt. ing(ir1,ir2)) nhbresmax=ing(ir1,ir2)
        rmsd2d(ir1,ir2)=rmsd2d(ir1,ir2)+nhbdist(ihb)
        if (ir1 .ne. ir2) then
          ing(ir2,ir1)=ing(ir1,ir2)
          rmsd2d(ir2,ir1)=rmsd2d(ir1,ir2)
        end if
      end do
      write (6,1017) bondname(1:lbondname),nhbresmax
      write (iout,1017) bondname(1:lbondname),nhbresmax
      write (iout,1010) bondname(1:lbondname)
      do jr=nreshb,1,-1
        write (iout,1003) jr,it3(jr),(ing(ir,jr),ir=1,nreshb)
      end do
      write (iout,1013) bondname(1:lbondname)
      do jr=nreshb,1,-1
        do ir=1,nreshb
          rmsd2d(ir,jr)=rmsd2d(ir,jr)/float(nframe)
          a1(ir)=rmsd2d(ir,jr)
        end do
        call writesmall(iout,a1,nreshb,jr,it3(jr))
      end do
c     Calculate average res-res bond percentage
      write (iout,1012) bondname(1:lbondname)
      write (iout,1014)(it3(ir),ir=1,nreshb)
      do jr=nreshb,1,-1
        call zeroit(a1,nreshb)
        do ir=1,nreshb
          if (ing(ir,jr) .gt. 0) a1(ir)=rmsd2d(ir,jr)/float(ing(ir,jr))
        end do
        call writesmall(iout,a1,nreshb,jr,it3(jr))
        if (ans .eq. 'a') call trnsfr(rmsd2d(1,jr),a1,nreshb)
      end do
c     Generate contact statistics ignoring multiple contacts
c     Sort contact list
      call sortbondlist(ixres,nbfound,indexbond,maxrsd,mxbonds)
      irescodeprev=0
      do ifr=1,nframe
        call readbitc(ires(1,ifr),ibond,nbfound,30,MAXITEMS)
        do ib=1,nbfound
          ihb=indexbond(ib)
          if (ibond(ihb) .gt. 0) then
            ir1=it2(ixres(ihbpair(1,ihb)))
            ir2=it2(ixres(ihbpair(2,ihb)))
            irescode=ir1*MAXRSD+ir2
            if (irescode .ne. irescodeprev) then
              irescodeprev=irescode
              ing(ir1,ir2)=ing(ir1,ir2)+1
              ing(ir2,ir1)=ing(ir1,ir2)
            end if
          end if
        end do
      end do
      write (iout,1016) bondname(1:lbondname)
      write (iout,1014)(it3(ir),ir=1,nreshb)
      do ir=1,nreshb
        do jr=1,nreshb
          ing(ir,jr)=0
        end do
      end do
      ip=1
      do while (ihb_pair_res(3,ip) .gt. 0)
        ir1=it2(ihb_pair_res(1,ip))
        ir2=it2(ihb_pair_res(2,ip))
        if (ir1*ir2 .gt. 0) ing(ir1,ir2)=ihb_pair_res(3,ip)
        ip=ip+1
      end do
      n_res_res=ip-1
      do jr=nreshb,1,-1
        call zeroit (a1,nreshb)
        do ir=1,nreshb
          a1(ir)=float(ing(ir,jr))/float(nframe)
        end do
        call writesmall(iout,a1,nreshb,jr,it3(jr))
        if (ans .eq. 'i') call trnsfr(rmsd2d(1,jr),a1,nreshb)
      end do
      rhbmax=0.0
      do ir=1,nreshb
        do jr=1,nreshb
          if (rhbmax .lt. rmsd2d(ir,jr)) rhbmax=rmsd2d(ir,jr)
        end do
      end do
      if (nreshb .gt. 0) then
        inc=max0(1,500/nreshb)
        iydel=150
c       call indexit(irrix,1,nreshb,0)
        call indexit(irrix,1,n_res_res,0)
        write (question,2002) bondname(1:lbondname)
        lq=38+lbondname
        call getreal(question,lq,rhbmax,rhbscalemax,0,74)
        call getint('Last residue to plot on the Y axis',34,
     -    nreshb,1,nreshb,iylst,94)
        call getint('First residue to plot on the Y axis',35,
     -    1,1,nreshb,iyfst,94)
        iyinc=iyfst-1
        nreshby=iylst-iyinc
        if (nreshby .le. 40) then
          nyclab=nreshby
          lyclab=1
c         Gather residue names and numbers
          do iy=iyfst,iylst
            ia=ifres(it3(iy))
            write (yclab(iy-iyinc),1008)
     -        line(index(ia))(irescol1:irescol2),iresno(ia)
            if (irescol2-irescol1+7 .gt. lyclab)
     -        lyclab=irescol2-irescol1+7
c           print *,iy,' YCLAB=',yclab(iy-iyinc)(1:lyclab)
          end do
        else
          nyclab=1
          lyclab=0
        end if
        scalefac=1.0
        iytop=0
        call plotmat(iplot,ing,rmsd2d,dc,nreshb,nreshby,0,0,0,iyinc,
     -    1,0,40,iydel,00,iytop,0.0,rhbscalemax,ncolcode,maxcolcode,
     -    ixdelsh,iydelsh,inc,scalefac,ixshuffle,ixshuffle,irrix,
     -    title,ltitle,'Residue-residue bond strength',29,1,' ',1,a1,
     -    yclab,nyclab,lyclab,1,MAX2D,1,MAX2D,MAX2D,ipspage,0)
c       Draw boundary between solute molecules
        write (iplot,*) '% Drawing separator lines'
        call rgbcolor(iplot,9)
        write (iplot,*) '2 lw'
        if (nreshby .gt. 1) scalefac=scalefac*0.95
        do i=1,nreshb-1
          if (isegno(ifres(it3(i+1))) .ne. isegno(ifres(it3(i)))) then
c           print *,'BOUND at i=',i,' ir=',it3(i),' ia=',ifres(it3(i+1))
            write (iplot,*) 'np'
            ixdell=ixdelsh
            if (nreshb .lt. 550) ixdell=ixdelsh+40
            write (iplot,1005) ixdell+i*inc,iydelsh
            write (iplot,1006) ixdell+i*inc,iydelsh+nreshb*inc
            write (iplot,1005) ixdell,iydelsh+(i-iyfst+1)*inc
            write (iplot,1006) ixdell+nreshb*inc,iydelsh+(i-iyfst+1)*inc
            write (iplot,*) 'sk'
          end if
        end do
        if (nreshb .ge. 5) iydel=iydel-50
        if (nreshb .lt. 5) iydel=iydel-20
        ixcent=max0(0,(nreshb*inc-50*ncolcode)/2)
        call colcodeminmax(iplot,25+ixcent,-iydel-5,0,ncolcode,
     -    maxcolcode,0.0,rhbscalemax)
        iydel=iydel-35
        write (iplot,1005) 75,iydel
        write (iplot,1004) bondname(1:lbondname),0.0,rhbscalemax
        iydel=iydel-15
        write (iplot,1005) 75,iydel
        write (iplot,1007) bondname(1:lbondname)
        iydel=iydel-15
        write (iplot,1005) 75,iydel
        if (ans .eq. 'i') then
          write (iplot,1009) 'Multiple contacts are ignored'
        else if (ans .eq. 'a') then
          write (iplot,1009) 'Multiple contacts are averaged'
        else
          write (iplot,1009) 'Multiple contacts are accumulated'
        end if
      else
        Print *,'No ',bondname(1:lbondname),'-bonding residue was found'
      end if
      write (iplot,*) 'showpage'
      print *,'END RES_RES_BOND'
      return
1000  format(' Number of different ',a,' bonds used=',i6,/,
     -  ' Number of frames analyzed=',i6,/)
1001  format(' List of the original indices of the ',a,
     -  '-bonding residues',/,(' Residue index(on plot)=',i4,
     -  ' residue #=',i5,' segment #=',i4,' residue index=',i5))
1003  format(' iy=',i4,' iy(orig)=',i5,10i5,/(23x,10i5))
1004  format('( Range of the ',a,'-bond fraction scale:',f8.3,'  -',
     -  f8.3,') show')
1005  format(2i5,' m')
1006  format(2i5,' l')
1007  format('( Residue-Residue ',a,' bond fractions ) show')
1008  format(1x,a,i5)
1009  format('( ',a,' ) show')
1010  format(/,' Number of different ',a,' bonds between residue pairs')
1012  format(/,' Fraction of time the ',a,' bonds are formed between ',
     -  'residue pairs',/,
     -  ' averaged over all bonds between these residues')
1013  format(/,' Average number of ',a,' bonds between residue pairs')
1014  format(' iy(orig):',13x,10i5,/,(23x,10i5))
c1015  format('( Upper-left triangle: iy is donor; ',
c     -  'Lower-right triangle: ix is donor ) show')
1016  format(/,' Fraction of time any ',a,' bond was formed between ',
     -  'residue pairs')
1017  format(/,' Maximum number of ',a,' bonds between two residues=',
     -  i4)
1018  format(' Number of residues involved in ',a,' bonding=',i5)
1019  format(' NOTE: residue numbers below refer to the residue ',
     -  'INDICES printed above')
1021  format(/,' === Residue-residue bond matrix data:')
2000  format('Do you want only ',a,'-bonded residues in the plot')
2002  format('Limit of the residue-residue ',a,' bond map scale')
      end
