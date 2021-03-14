      subroutine unpacktraj(inptrajtyp,mmctrajtyp,inpcrdtyp,iotyp,
     -  inpcrdtyporg,n,nslt,naslv,islvw,numres,numslv,resnamslv,iatnum,
     -  ifchrg,ixres,atnames,resnames,segnames,altcol,inscol,ninsres,
     -  marker,ireseq,iresnrestart,iresidrestart,c,rprox,cv,charge,
     -  ntitlin,title,nneig,nneiga,nhbneig,ineig,nhneig,nnneig,ncneig,
     -  nsneig,npneig,isegno,iresno,line,blankline,ifree,icntrl,innlist,
     -  mmtype,ibnd,index,indexn,indexo,molresflag,hblimfac,angmin,
     -  ioppbc,cell,ncell,ionefile,imodel,numsel,iconfsel,idupl,
     -  itemp1,nclstmem,numconsec,incr_fileno,nextconfsel,iqspaceask,
     -  keeprem,radtodeg,ioverallconf,lentest,maxconfsel,maxrepconf,
     -  maxng,maxrsd,maxrec,mx2d)
      character* 132 line(maxrec),blankline
      dimension nneig(maxrec),ineig(maxng,maxrec),iatnum(maxrec),
     -  ifchrg(maxrec),ixres(maxrec),c(3,maxrec),rprox(maxrec),
     -  cv(maxrec),nhbneig(maxrec),nneiga(maxrec),isegno(maxrec),
     -  iresno(maxrec),ifree(maxrec),charge(maxrec),nhneig(maxrec),
     -  nnneig(maxrec),ncneig(maxrec),nsneig(maxrec),npneig(maxrec),
     -  mmtype(maxrec),ibnd(maxng,maxrec),index(maxrec),indexn(maxrec),
     -  indexo(maxrec),molresflag(maxrsd),iconfsel(maxconfsel),
     -  idupl(maxrsd),itemp1(maxrec),nclstmem(mx2d),cell(3,ncell)
      real*8 etotsum,etotsum2
      character*1 altcol(maxrec),inscol(maxrec)
      character*4 pflsv(100),segnames(maxrsd)
      character*8 resnamslv,atnames(maxrec),resnames(maxrsd),namesv(100)
      character*6 marker(16)
      character*80 trtitle(32),title
      character*200 inpfile
      dimension icntrl(20)
      common /logging/ logfile,ipredict
      character*5 crdext
      common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
     -  iommod,iommc,iommc4,iogro,iomol2,iomae,iocif,ioxxx,ioins,
     -  ionxyz,iosxyz,iosxyzrq,iograsp,iofull,lext(19),crdext(19)
      character*200 trajnam,outfile
      common /trajname/ trajnam,outfile,ltrajnam,namleno,ifirsttraj,
     -  ifirsttraj2,ilasttraj,ilasttraj2,incrementtraj,incrementtraj2
      character*11 trajformatname
      common /trajectory/ nmmccheck,iftrajtyp(6),trajformatname(6)
      dimension iasv(100),qsv(100)
      data ifirst /1/,ilast /999999/,increment /1/,inpt /70/
c     print *,'UNPACKTRAJ inptrajtyp,mmctrajtyp=',inptrajtyp,mmctrajtyp
      nmax0=n
      nmax=n
      etotsum=0.d0
      etotsum2=0.d0
      etotmin=1.e+30
      etotmax=-etotmin
      emax=0.0
      netot=0
      neskip=0
      neavskip=0
      ielim=1
      if (numsel .gt. 0) write (6,1010) (iconfsel(i),i=1,numsel)
      call opentraj(c,1,inpt,inptrajtyp,n,ntitltr,trtitle,
     -  inpcrdtyp,ifirst,ilast,increment,maxconf,ninconf,noutconf,
     -  natom,nfreeat,ifree,icntrl,1,mmctrajtyp,trajnam,ltrajnam,
     -  'trajectory',10,iconfsel,numsel,1,0,0,0,icellfound,notfnd,0,
     -  0,lentest,0,0,0,maxconfsel,maxrec)
      if (ionefile .eq. 0) then
        call getname(inpfile,namleni,
     -    'Output coordinate file (without the extension)',46,200,'',0,
     -    0,0,0)
        ic=namleni+1
        inpfile(ic:ic)='.'
        inpfile(ic+1:ic+3)=crdext(iotyp)
        namleni=namleni+4
      else
        call getname(inpfile,namleni,'Output coordinate file',22,200,
     -    '',0,0,0,0)
      end if
      ifail=0
      iaskatnum=0
      nntest=0
      nclustermem=0
      lastwrite=0
c     write (6,7866) numsel,(iconfsel(i),i=1,numsel)
c7866 format(' NUMSEL=',i4,' ICONFSEL=',20i6)
      do while (ninconf .lt. ilast .and. lastwrite .eq. 0 .and.
     -          ifail .eq. 0)
c       Read a conformation
        call readtraj(inpt,inptrajtyp,mmctrajtyp,n,naslv,nwatr,nslt,
     -    trajformatname(inptrajtyp),ntitltr,trtitle,trajnam,ltrajnam,
     -    natom,nfreeat,ifree,icntrl,c,ninconf,noutconf,
     -    increment,inpcrdtyp,ietot,etot,ifail,ifirst,ilast,iconfsel,
     -    numsel,maxrepconf,nmc,lentest,1.0,maxconf,maxconfsel,maxrec)
        if (ifail .eq. 0) then
          if (nntest .eq. 0) then
            call comparetop(c,n,nneig,ineig,iatnum,innlist,n,1,
     -        cell,ncell,ioppbc,maxng,maxrec)
            nntest=1
          end if
          ifilenum=0
          if (numsel .gt. 0) then
c           List-directed unpacking - check list
            if (ninconf .eq. iconfsel(nextconfsel)) then
c             Selection found - increment nextconfsel
              ifilenum=ninconf
              nclustermem=nclstmem(nextconfsel)
              if (nextconfsel .eq. numsel) lastwrite=1
              if (nextconfsel .lt. min0(numsel,maxconfsel))
     -          nextconfsel=nextconfsel+1
            else if (iconfsel(nextconfsel) .eq. 0) then
              ifail=1
            end if
          else
c           Full unpacking
            if (ninconf .ge. ifirst .and.
     -        mod(ninconf-ifirst,increment) .eq. 0)
     -          ifilenum=noutconf+incr_fileno+1
          end if
          if (ietot .eq. 1) then
            if (ielim .eq. 1) then
              call askyn(
     -          'Do you want to give an energy threshold (max)',45,1,-1,
     -          iemax,0,0)
              if (iemax .eq. 1)
     -          call getreal('Energy maximum to use',21,0.0,emax,0,0)
              ielim=0
            end if
            if (iemax .eq. 1 .and. etot .gt. emax) then
              if (ifilenum .gt. 0) neskip=neskip+1
              ifilenum=0
            end if
          end if
          if (ifilenum .gt. 0) then
c           Write a conformation
            noutconf=noutconf+1
            ifnumw=ifilenum
            if (numconsec .gt. 0) ifnumw=noutconf
            if (noutconf .eq. 1 .or. ionefile .eq. 0) then
              if (ionefile .eq. 0) then
                call filenamenum(inpfile,namleni,outfile,nl1,ifnumw,+2)
              else
                outfile=inpfile
                nl1=namleni
              end if
              call openfile(20,0,'output',6,'new',outfile,nl1,
     -          notfnd,2,1,1,0,ioverallconf)
              if (notfnd .eq. 1) write (6,1012) outfile(1:nl1)
              if (notfnd .eq. 1 .or. 
     -            (ipredict .eq. 1 .and. ioverallconf .eq. 0)) then
                call askyn(
     -            'Do you want to overwrite all existing files',43,1,-1,
     -            ioverallconf,0,0)
                if (ioverallconf .eq. 0) stop
                call openfile(20,0,'output',6,'new',outfile,nl1,
     -            notfnd,2,1,1,0,ioverallconf)
              end if
            end if
            if (n .gt. nmax) then
c             Add extra solvents
              nmax0=nmax
              call updatesolvents(iaskatnum,1,nmax,n,naslv,iasv,namesv,
     -          resnamslv,qsv,pflsv,iwat,iatnum,noutconf,iotyp,maxrec)
            end if
            call writeconf(20,inpcrdtyp,iotyp,inpcrdtyporg,nmax0,n,nslt,
     -        naslv,islvw,iasv,namesv,qsv,pflsv,1,1,imodel,iatnum,
     -        ifchrg,noutconf,innlist,c,rprox,cv,ixres,iresno,atnames,
     -        resnames,segnames,charge,isegno,altcol,inscol,ninsres,
     -        marker,ntitlin,0,title,ireseq,iresnrestart,iresidrestart,
     -        nneig,nneiga,nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,
     -        npneig,numres,numslv,resnamslv,line,blankline,mmtype,ibnd,
     -        index,indexn,indexo,0,molresflag,idupl,itemp1,hblimfac,
     -        angmin,0,1,0,0,noutconf-1,5,iqspaceask,0,ifilenum,etot,
     -        ietot,nclustermem,0,keeprem,iwriteatsym,radtodeg,
     -        maxrepconf,maxng,maxrsd,maxrec)
            if (ietot .eq. 1) then
              if (etot .lt. emax) then
                etotsum=etotsum+etot
                etotsum2=etotsum2+etot**2
                if (etot .lt. etotmin) etotmin=etot
                if (etot .gt. etotmax) etotmax=etot
                netot=netot+1
              else
                neavskip=neavskip+1
              end if
            end if
            if (ionefile .eq. 0) close (20)
            if (noutconf .eq. 1) inpcrdtyp=iotyp
          end if
        end if
      end do
      if (netot .gt. 0) then
        eav=etotsum/dfloat(netot)
        sd=dsqrt(etotsum2/dfloat(netot)-(etotsum/dfloat(netot))**2)
        write (6,1021) eav,sd,etotmin,etotmax,netot
        if (neavskip .gt. 0 .and. iemax .eq. 0) write (6,1022) neavskip
      end if
      if (ionefile .eq. 1) close (20)
      close(inpt)
      write (6,1011) noutconf
      if (neskip .gt. 0) write (6,2127) neskip,emax
      return
1010  format(' Unpacking trajectory frames ',(/,10i7))
1011  format(' Completed unpacking trajectory into ',i6,
     -  ' configurations')
1012  format(' Problem opening file ',a)
1021  format(' <E>=',f10.4,' S.D.=',f10.4,' Range: [',f10.4,',',
     -  f10.4,'] (N=',i5,')')
1022  format(' Note:',i5,' structures with non-negative energy ',
     -  'were skipped for the statistics')
2127  format(' Note:',i5,' structures with energy >',e12.5,
     -  'were skipped')
      end
