      subroutine maketraj(nconfig,maxconf,c,c1,rprox,cv,ih,nslt,n,naslv,
     -  islvw,ntitltr,trtitle,line,altcol,inscol,ninsres,inpcrdtyp,
     -  newinp,inptrajtyp,ioutrajtyp,mmctrajtyp,iout,iatnum,ifchrg,
     -  innlist,iresno,ixres,atnames,resnames,segnames,charge,isegno,
     -  marker,ntitlin,title,ireseq,iresnrestart,iresidrestart,nneig,
     -  nneiga,nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,npneig,numres,
     -  numslv,resnamslv,blankline,mmtype,ibnd,index,indexn,indexo,
     -  indexs,icntrl,icntrlr,ifree,molresflag,idupl,itemp1,ireorder,
     -  icellfound,natsdel,nwrmax,hblimfac,angmin,iqspaceask,keeprem,
     -  iwriteatsym,radtodeg,etot,noheadchange,maxrepconf,maxng,mxres,
     -  maxrec)
      dimension icntrl(20),icntrlr(20),ifree(n)
      character*1 altcol(maxrec),inscol(maxrec)
      character*80 trtitle(32)
      character* 132 line(maxrec),blankline
      character*80 title
      character*4 segnames(mxres)
      character*8 resnamslv,atnames(maxrec),resnames(mxres)
      character*6 marker(16)
      dimension c(3,maxrec),c1(3,maxrec),rprox(maxrec),cv(maxrec),
     -  ih(maxrec),nneig(maxrec),ineig(maxng,maxrec),iatnum(maxrec),
     -  ifchrg(maxrec),charge(maxrec),nhbneig(maxrec),nneiga(maxrec),
     -  nhneig(maxrec),nnneig(maxrec),ncneig(maxrec),nsneig(maxrec),
     -  npneig(maxrec),isegno(maxrec),iresno(maxrec),ixres(maxrec),
     -  mmtype(maxrec),ibnd(maxng,maxrec),index(maxrec),indexn(maxrec),
     -  indexo(maxrec),indexs(maxrec),molresflag(mxres),idupl(mxres),
     -  itemp1(maxrec)
      character*5 crdext
      common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
     -  iommod,iommc,iommc4,iogro,iomol2,iomae,iocif,ioxxx,ioins,
     -  ionxyz,iosxyz,iosxyzrq,iograsp,iofull,lext(19),crdext(19)
      real*8 xtlabc,xtlabc0
      common /boxdat/ xtlabc(6),xtlabc0(6),box(3),box0(3),edge_gen(3,3),
     -  cell0(3,27),cell(3,27),cellalt(3,27),
     -  ncell,ioppbc,noboxinfoar,noboxinfow,noboxrep,
     -  istuner,iboxtypfound,ixcrd(3),ixang(3),ixyzhex(3),
     -  ixyzhextraj(3),isizewarn
      character*1 xyz
      common /axislab/ xyz(3)
      common /saveinfo/ nwrite,nfreeatwrite
      character*4 pflsv(100)
      character*8 namesv(100)
      dimension cmin(3),cmax(3),shiftmmc(3),iasv(100),qsv(100)
      real*8 detot
      data ifsavdef /1/,incrtdef /1/,ifsavvdef/0/,ichvers /27/,
     -  nfixatdef /0/,tsfsdef /1.0/,ndegfree /0/
      character*4 charmmheader(2)
      character*11 charmmheadname(2)
      data charmmheader/'CORD','VELD'/
      data charmmheadname/'coordinates','velocities '/
c     print *,'MAKETRAJ nslt,n,naslv=',nslt,n,naslv
c     print *,'MAKETRAJ ioutrajtyp,ioppbc,icntrlr(11)=',
c    -  ioutrajtyp,ioppbc,icntrlr(11)
      limic11=1
      if (nconfig .eq. 1) then
        limic11=0
        nwrite=n-natsdel
        if (natsdel .gt. 0)
     -    print *,'Number of atoms written in a configuration=',nwrite
c       Create and write header
        if (ioppbc .eq. -1) then
          do k=1,3
            box0(k)=999.0
            box(k)=999.0
          end do
          noboxinfow=1
        else
          noboxinfow=0
        end if
c       Decide about cell info writing
        iboxdef=-1
        if (iboxtypfound .gt. 0) iboxdef=1
        if ((ioutrajtyp .eq. 1 .or. ioutrajtyp .eq. 3) .and.
     -       ioppbc .ge. 0)
     -    call askyn('Do you want box information written',35,
     -      0,iboxdef,noboxinfow,0,0)
        if (ioutrajtyp .eq. 1) then
c         Charmm
          if (noheadchange .eq. 0) then
            if (noboxinfow .eq. 0) then
              icntrl(11)=1
              if (ioppbc .ge. 0) then
                write (6,1104)
                call askyn('Do you want repeated box information',36,
     -            0,1,noboxrep,0,0)
                if (noboxrep .eq. 0) icntrl(11)=2
              end if
              if (noboxinfoar .eq. -1) then
c               Input box size
                print *,'Only rectangular box input is implemented'
                do k=1,3
                  call getreal(
     -              'Box dimension in the '//xyz(k)//' direction',32,
     -              999999.0,xk,1,0)
                  xtlabc(ixcrd(k))=xk
                  xtlabc(ixang(k))=90.0d0
                end do
              end if
            else
              icntrl(11)=0
              call zeroitd(xtlabc,6)
            end if
          else
            icntrl(11)=icntrlr(11)
          end if
        else if (ioutrajtyp .eq. 2) then
c         Amber
          noboxrep=1
          irepdef=-1
          if (iboxtypfound .gt. 1) irepdef=1
          if (ioppbc .ge. 0) then
            call askyn('Do you want repeated box information',36,
     -        0,irepdef,noboxrep,127,9)
          end if
          if (inptrajtyp .eq. 1. and. icntrl(11) .ne. -1) then
c           Use Charmm box size
            do k=1,3
              box(k)=xtlabc(ixcrd(k))
            end do
          else if (icellfound .eq. 0 .and. ioppbc .ge. 0) then
c           Input box size
            call getxyz('Initial box dimension in the ',29,
     -        ' direction',10,999999.0,box,1,0)
          end if
        end if
        if (ioutrajtyp .eq. 1) then
c         Charmm
          icntrl(1)=maxconf
          if (inptrajtyp .eq. 1) then
            incrtdef=icntrlr(2)
            ifsavdef=icntrlr(3)
            ifsavvdef=icntrlr(5)
            ndegfree=icntrlr(8)
            nfixatdef=icntrlr(9)
            ichvers=icntrlr(20)
c           call copybits(icntrlr(10),tsakma)
            call copyintgreal(icntrlr(10),tsakma,1)
            tsfsdef=tsakma*48.88821
          end if
          if (noheadchange .eq. 0) then
            call getint('Number of fixed atoms to put in the header',
     -        42,nfixatdef,1,n,icntrl(9),50)
            call getint(
     -        'Coordinate saving frequency to put in the header',48,
     -        ifsavdef,1,0,icntrl(3),50)
            call getint('Previous run steps to put in the header',39,
     -        incrtdef,1,0,icntrl(2),50)
            call getint(
     -        'Velocity saving frequency to put in the header',46,
     -        ifsavvdef,1,0,icntrl(5),50)
            call getint('Charmm version to put in the header',35,
     -        ichvers,1,40,icntrl(20),50)
            call getreal('Time step/fs to put in the header',33,tsfsdef,
     -        tsfs,1,0)
            tsakma=tsfs/48.88821
c           call copybits(tsakma,icntrl(10))
            call copyintgreal(icntrl(10),tsakma,0)
            icntrl(4)=maxconf*ifsavdef
            print *,'Number of creation steps=',icntrl(4)
            call askyn(
     -       'Do you want to create a velocity trajectory (VELD)',50,
     -       1,-1,icharmmheader,0,0)
            icharmmheader=icharmmheader+1
            write (6,1103) charmmheadname(icharmmheader),
     -        charmmheader(icharmmheader)
          else
            icntrl(2)=icntrlr(2)
            icntrl(3)=icntrlr(3)
            icntrl(5)=icntrlr(5)
            icntrl(8)=icntrlr(8)
            icntrl(9)=icntrlr(9)
            icntrl(10)=icntrlr(10)
            icntrl(20)=icntrlr(20)
            icharmmheader=1
          end if
          nfreeat=n-icntrl(9)
          nfreeatwrite=nfreeat
          nfadel=0
          if (natsdel .gt. 0) then
c           See if nfreeat has to be reduced when solvents are dropped
            if (ioutrajtyp .eq. 1 .and. icntrl(9) .gt. 0) then
              nfadel=0
              do ia=1,nfreeat
                if (indexs(ifree(ia)) .eq. 0) then
                  nfadel=nfadel+1
                else
                  ifree(ia-nfadel)=ifree(ia)
                end if
              end do
              nfreeatwrite=nfreeat-nfadel
            end if
          end if
          if (inptrajtyp .eq. 1) then
            if (natsdel .gt. 0) ndegfree=ndegfree-3*nfadel
          else
            ndegfree=nfreeatwrite*3-6
          end if
          icntrl(8)=ndegfree
          write (iout) charmmheader(icharmmheader),icntrl
          write (iout) ntitltr,(trtitle(i),i=1,ntitltr)
          write (iout) nwrite
          if (icntrl(9) .gt. 0) write (iout) (ifree(i),i=1,nfreeatwrite)
        else if (ioutrajtyp .eq. 2) then
c         Amber
          write (iout,1000) trtitle(1)
        else if (ioutrajtyp .eq. 3) then
c         MMC
c         if (naslv .eq. 3) then
c           call getint('Number of solvent atoms per solvent molecule',
c    -        44,999999,1,0,naslv,0)
c         end if
          call extension(c,ih,0,1,nslt,cmin,cmax,shiftmmc,0,0,v)
          write (6,1101) shiftmmc
        else if (ioutrajtyp .gt. 5) then
          print *,'PROGRAM ERROR: invalid trajectory type:',
     -      ioutrajtyp
        end if
      end if
      if (ireorder .gt. 0) then
c       Create array with new order
        ndel=0
        do i=1,nwrmax
          if (indexs(i) .gt. 0) then
            do k=1,3
              c1(k,i-ndel)=c(k,indexs(i))
            end do
          else
            ndel=ndel+1
          end if
        end do
        nwrite=nwrmax-ndel
      end if
      if (ioutrajtyp .eq. 1) then
c       Charmm
c       print *,'MAKETRAJ icntrl(11),limic11=',icntrl(11),limic11
        if (icntrl(11) .gt. limic11) write (iout) xtlabc
        if (ireorder .eq. 0) then
          if (icntrl(9) .eq. 0 .or. nconfig .eq. 1) then
            do k=1,3
              write (iout) (c(k,i),i=1,nwrite)
            end do
          else
            do k=1,3
              write (iout) (c(k,ifree(i)),i=1,nfreeatwrite)
            end do
          end if
        else
          if (icntrl(9) .eq. 0 .or. nconfig .eq. 1) then
            do k=1,3
              write (iout) (c1(k,i),i=1,nwrite)
            end do
          else
            do k=1,3
              write (iout) (c1(k,ifree(i)),i=1,nfreeatwrite)
            end do
          end if
        end if
      else if (ioutrajtyp .eq. 2) then
c       Amber
        if (ireorder .eq. 0) then
          write (iout,1001) ((c(k,i),k=1,3),i=1,nwrite)
        else
          write (iout,1001) ((c1(k,i),k=1,3),i=1,nwrite)
        end if
        if ((nconfig .eq. 1 .and. ioppbc .gt. 0) .or. noboxrep .eq. 0)
     -    write (iout,1001) box
      else if (ioutrajtyp .eq. 3) then
c       MMC
        if (mmctrajtyp .eq. 1) then
          detot=etot
          nwat=(n-nslt)/naslv
          write (iout) nwat,n,1.0d0,nconfig,0,
     -      0,0,1,nslt,detot,0.d0,0.d0,0.d0,0.d0,0.0
          if (ireorder .eq. 0) then
            write (iout) ((c(k,i),k=1,3),i=1,n)
          else
            write (iout) ((c1(k,i),k=1,3),i=1,n)
          end if
          if (noboxinfow .eq. 0) write (iout) box
        else
          print *,'Non-binary MMC trajectory write is not implemented'
          stop
        end if
      else if (ioutrajtyp .eq. 4 .or.
     -         (ioutrajtyp .eq. 5 .and. nconfig .eq. 1)) then
c       Macromodel
        inptp=inpcrdtyp
        inpcrdtyporg=inpcrdtyp
        if (nconfig .gt. 1 .and. newinp .eq. 0) then
          inptp=iommod
        end if
        call writeconf(iout,inptp,iommod,inpcrdtyporg,nwrite,nwrite,
     -    nslt,naslv,islvw,iasv,namesv,qsv,pflsv,1,1,0,iatnum,ifchrg,
     -    nconfig,innlist,c,rprox,cv,ixres,iresno,atnames,resnames,
     -    segnames,charge,isegno,altcol,inscol,ninsres,marker,ntitlin,0,
     -    title,ireseq,iresnrestart,iresidrestart,nneig,nneiga,nhbneig,
     -    ineig,nhneig,nnneig,ncneig,nsneig,npneig,numres,numslv,
     -    resnamslv,line,blankline,mmtype,ibnd,index,indexn,indexo,0,
     -    molresflag,idupl,itemp1,hblimfac,angmin,0,1,0,0,nconfig-1,5,
     -    iqspaceask,0,1,0.0,0,0,0,keeprem,iwriteatsym,radtodeg,
     -    maxrepconf,maxng,mxres,maxrec)
      else
c       Just write the shortened version (Xcluster)
        write (iout,1010) -n,title,((c(k,i),k=1,3),i=1,n)
      end if
      return
1000  format(a80)
1001  format(10f8.3)
1010  format(i5,a,/,(5x,3f12.5))
1101  format(' The system will be shifted by ',3f10.5)
1103  format(' Trajectory is assumed to contain ',a,' (header:',a,')')
1104  format(' For Charmm 27 and higher versions cell information can ',
     -  'be stored ',/,' for each frame. This is required for (T,P,N) ',
     -  'runs')
      end
