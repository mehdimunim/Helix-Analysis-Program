      subroutine readconf(inpt,inpcrdtyp,iruntyp,iruntypclean,
     -  inpfile,namleni,outfile,namleno,line,title,ititle,trtitle,
     -  ntitlin,ntitltr,nconfig,ncol,n,nlines,nsegm,nosegid,iisegcol,
     -  index,c,iresno,iatnum,isegno,froccin,charge,segid4,altcol,
     -  inscol,iundef,naltnam,nrecdel,ninsres,lineread,nneig,ineig,iha,
     -  itrunc,iendfound,ioktoend,istopatend,iqha,imodelkeep,newsega,
     -  newsegr,neednnlist,ietotread,etotread,bfacavg,molname,
     -  molnamelen,pdbid,iclone,iask77_78,iread77_78,iwarn77_78,
     -  iblank67_76,iask67_76,innlistread,iout,maxrepconf,maxng,maxrec,
     -  maxrsd)
      character*(*) inpfile,outfile
      character* 132 line(maxrec)
      character*80 title,trtitle(32),molname,linewr
      character*4 segid4(maxrsd),pdbid
      character*1 altcol(maxrec),inscol(maxrec)
      dimension c(3,maxrec),index(maxrec),iresno(maxrec),iatnum(maxrec),
     -  isegno(maxrec),froccin(maxrec),charge(maxrec),iisegcol(2,13),
     -  nneig(maxrec),ineig(maxng,maxrec),bfacavg(maxrsd)
      common /clonedat/ nclone,iaclnf(1000),iaclnl(1000),ncopcln(1000)
      character*11 formatname
      common /formats/ iqdconv(20),formatname(19)
      character*5 crdext
      common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
     -  iommod,iommc,iommc4,iogro,iomol2,iomae,iocif,ioxxx,ioins,
     -  ionxyz,iosxyz,iosxyzrq,iograsp,iofull,lext(19),crdext(19)
      character*2 iatnm2
      common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
     -  mmatno(64),iatnm2(99)
      common /askCA/ iaskca,nocalcium
      common /logging/ logfile,ipredict
      real*8 bfacsum
      character*1 altnam(50),s1,asterisk,ans
      character*2 atsymbol
      character*4 s4,segidprev
      character*8 atomnam
      character*8 segmid
      character*6 rnu
      character*51 question
      character* 132 lineinp,lineprev,linenext
      character*200 altfile
      dimension noblank3(132),in12(2)
      data nnoblank3 /15/, noblank3
     -  /1,5,11,13,19,21,27,29,35,37,43,45,51,119,124,117*0/
      data asterisk/'*'/, question
     -   /'Clone   : original index of the first and last atom'/
      iaskagain=0
      nclone=0
      iendfound=0
      iunrecog=0
      naltnam=0
      itrunc=0
      icinc=0
      nosegid=0
      etotread=0.0
      bfacsum=0.d0
      nbfacsum=0
      nbfacsumtot=0
      ikeepfullalt=-1
      innlistread=0
      naltrec=0
      naltdel=0
      ninsres=0
      call setcol(inpcrdtyp,ncol,idcol,ialtcol,iinscol,
     -  inamcol1,inamcol2,irescol1,irescol2,iccol1,iccol2,
     -  iresncol1,iresncol2,iseqncol1,iseqncol2,isegcol1,isegcol2,
     -  iresidcol1,iresidcol2,iqcol1,iqcol2,ipotcol1,ipotcol2,
     -  iocccol1,iocccol2,ichemcol1,ichemcol2,nrescol,nresncol,
     -  nsegcol,nnamcol,iofull)
      iresnoprev=0
      iresinc=0
      incresno=0
      lineread=0
      linefread=0
      if (iresncol2-iresncol1 .lt. 8)
     -  incresno=10**(iresncol2-iresncol1+1)-1
      nsegmid=isegcol2-isegcol1+1
100   if (ischarmm(inpcrdtyp) .eq. 1) then
c-------Charmm input
        call blankout(lineinp,1,ncol)
        do i=1,maxrec
          if (linefread .eq. 0) read (inpt,1001,end=9999) lineinp
          lineread=0
          line(i)=lineinp
          if (i .eq. 1) then
             if (lineinp(1:1) .ne. '*') then
               write (6,1215) inpfile(1:namleni),'Charmm CRD'
               namleni=0
               iunrecog=1
             end if
             title=lineinp(3:80)
          end if
          if (lineinp(1:1) .ne. '*') then
            go to 11
          else if (i .gt. 32) then
            write (6,1213) i
            stop
          else
            trtitle(i)=lineinp(3:80)
            call checkforetot(1,lineinp,nconfig,etotread,ietotread,1)
          end if
          ntitlin=i
          ntitltr=i
          if (nconfig .le. maxrepconf) then
            call lastchar(lineinp,ifc,ncol)
            write (6,1010) lineinp(1:ifc)
          end if
        end do
c       Decode number of atoms
11      lineread=ntitlin+1
        call lastchar(lineinp,lc,ncol)
        if (lc .le. 1 .and. lineinp(1:1) .eq. ' ') go to 9917
        if (inpcrdtyp .eq. iocha) read (lineinp(1:6),*,err=9915) n
        if (inpcrdtyp .eq. iochaex) read (lineinp(1:10),*,err=9915) n
        line(lineread)=lineinp
        nlines=ntitlin+1+n
        if (n .lt. 1) then
          print *,'ERROR: non-positive number of atoms found:',n
          stop
        end if
        call checkdim(nlines,MAXREC,'MAXREC',6,'number of lines',15,0)
        print *,'Reading',n,' atoms; configuration #:',nconfig
        segmid='        '
        do i=1,n
          read (inpt,1001,end=3000) lineinp
          line(ntitlin+1+i)=lineinp
          index(i)=ntitlin+1+i
          if (i .gt. 1) then
            if (lineinp(isegcol1+icinc:isegcol2+icinc) .ne.
     -          segmid(1:nsegmid)) then
              nsegm=nsegm+1
              iresinc=0
            end if
          end if
          isegno(i)=nsegm
          call readint(lineinp,iresncol1+icinc,iresncol2+icinc,
     -      iresno(i),2,0,ierr)
          if (ierr .gt. 0) go to 9919
          if (iresnoprev-iresno(i) .eq. incresno .and. incresno .ne. 0)
     -      iresinc=iresinc+incresno+1
          iresnoprev=iresno(i)
          iresno(i)=iresno(i)+iresinc
          lineread=ntitlin+i+1
          if (inpcrdtyp .eq. iocha) then
            read (lineinp(iccol1+icinc:iccol2+icinc),1101,err=9911)
     -        (c(k,i),k=1,3)
          else
            read (lineinp(iccol1+icinc:iccol2+icinc),1105,err=9911)
     -        (c(k,i),k=1,3)
          end if
          if (iqcol2 .ge. iqcol1)
     -      read (lineinp(iqcol1:iqcol2),*,err=9912,end=8801) charge(n)
          go to 8802
c         Blank charge field should be zero
8801      charge(n)=0.0
8802      iatnum(i)=
     -      ianum(lineinp(inamcol1+icinc:inamcol2+icinc),1,nnamcol)
          lineprev=lineinp
          if (inpcrdtyp .eq. iocha) then
            segmid(1:nsegmid)=lineinp(isegcol1+icinc:isegcol2+icinc)
          else
            call leftadjustline(lineinp,isegcol1,isegcol2)
            segmid(1:4)=lineinp(isegcol1:isegcol1+3)
          end if
          if (i .eq. 999999) icinc=1
        end do
      else if (ispdb(inpcrdtyp) .gt. 0) then
c-------PDB input
        iundef=0
        ititfound=0
        ititrfound=0
        ntitltr=0
        ntitlin=0
        linenext_read=0
        namlena=0
        molnamelen=0
        iatseqnumprev=0
        ignoresegnojump=0
        if (inpcrdtyp .ne. iobpdb) then
           iask67_76=0
           iblank67_76=0
        end if
        lineinp(1:1)='X'
        if (nconfig .eq. 0 .and. inpcrdtyp .ne. ioa4pdb) then
          if (ipredict .eq. 1) then
            write (6,2077)
            call askyn(
     -        'Do you want to read chemical symbols from col 77-78',51,
     -        1,1,iread77_78,0,0)
            if (iread77_78 .eq. 1) nocalcium=0
            iask77_78=0
          end if
        else
          iask77_78=0
        end if
        i=0
        do while (lineinp(1:4) .ne. 'ATOM' .and.
     -            lineinp(1:6) .ne. 'HETATM' .and. i .lt. maxrec)
          i=i+1
          call blankout(lineinp,1,ncol)
          read (inpt,1001,end=9999) lineinp
          call lastchar(lineinp,lcinp,ncol)
          lineread=i
          call checkforetot(6,lineinp,nconfig,etotread,ietotread,1)
          line(i)=lineinp
          if (lineinp(1:5) .eq. 'TITLE' .and. ititrfound .eq. 0) then
            title=lineinp(7:80)
            ititle=1
            ititrfound=1
            ititfound=1
          end if
          if (lineinp(1:6) .eq. 'HEADER') pdbid=lineinp(63:66)
          if (ititfound .eq. 0 .and. (lineinp(1:6) .eq. 'HEADER'
     -        .or. lineinp(1:6) .eq. 'REMARK')) then
            title=line(i)(7:80)
            ititfound=1
          end if
          if ((lineinp(1:5) .eq. 'TITLE' .or. lineinp(1:6) .eq. 'HEADER'
     -        .or. lineinp(1:6) .eq. 'REMARK')
     -        .and. ntitltr .lt. 32) then
            ntitltr=ntitltr+1
            trtitle(ntitltr)=lineinp(7:80)
          end if
          if (nconfig .le. maxrepconf) then
            if (lineinp(1:6) .eq. 'HEADER' .or.
     -          lineinp(1:5) .eq. 'TITLE' .or.
     -          lineinp(1:6) .eq. 'REMARK' .and. i .lt. 11) then
              write (iout,1015) lineinp(1:min0(79,lcinp))
            end if
          end if
          if (lineinp(1:31) .eq. 'REMARK    Name               : ') then
            ic=32
            call nextblank(lineinp,ic,ncol)
            molnamelen=ic-32
            molname(1:molnamelen)=lineinp(32:ic-1)
          end if
          if (lineinp(1:5) .eq. 'MODEL' .or. ipredict. eq. 1) then
            if (imodelkeep .eq. -1) call quiz(ans,imodelkeep,'k',' ',0,
     -        'MODEL record treatment',22,0,5,6,0)
            if (imodelkeep .gt. 1) then
              line(lineread)(1:6)='REMARK'
            end if
          end if
          ntitlin=ntitlin+1
        end do
        index(1)=ntitlin
        nlines=ntitlin
        ntitlin=ntitlin-1
        n=0
        nhdel=0
        nudel=0
        lineread=ntitlin
        do while (.true.)
          if (lineread .gt. ntitlin) then
            if (linenext_read .eq. 0) then
              call blankout(lineinp,1,ncol)
              read (inpt,1001,end=23) lineinp
            else
              lineinp=linenext
              linenext_read=0
            end if
            call checkdim(nlines,MAXREC,'MAXREC',6,'number of lines',15,
     -        0)
            line(lineread+1)=lineinp
            nlines=nlines+1
          end if
          lineread=lineread+1
          if (lineinp(1:5) .eq. 'MODEL') then
            if (imodelkeep .gt. 1) line(lineread)(1:6)='REMARK'
          else if (lineinp(1:6) .eq. 'ENDMDL') then
            if (imodelkeep .eq. 2) then
              line(lineread)(1:6)='REMARK'
            else if (imodelkeep .eq. 3) then
              line(lineread)(1:6)='END   '
            else if (imodelkeep .eq. 4) then
              line(lineread)(1:6)='TER   '
            end if
            go to 23
          else if (lineinp(1:3) .eq. 'END' .and. lineinp(4:7).ne. 'ROOT'
     -      .and. lineinp(4:9) .ne. 'BRANCH') then
c           Don't check for data after END for trajectory combining
            if (iruntyp .eq. 12 .and. lineinp(1:3) .eq. 'END') go to 23
            if (istopatend .eq. 1) then
              go to 23
            else if (istopatend .eq. -1) then
c             See if file ends
              call blankout(linenext,1,ncol)
              read (inpt,1001,end=23) linenext
              call lastchar(linenext,ilc,ncol)
              print *,'Data after END was found:'
              print *,linenext(1:ilc)
              linenext_read=1
              call askyn(
     -          'Do you want to replace middle END records with TER',
     -          50,0,-1,istopatend,000,0)
              if (istopatend .eq. 1) go to 23
              lineinp(1:3)='TER'
            else
c             Replace END with TER
              lineinp(1:3)='TER'
            end if
          end if
          if (lineinp(1:3) .eq. 'TER' .and. nbfacsum .gt. 0) then
            if (nbfacsum .gt. 0) bfacavg(nsegm)=bfacsum/nbfacsum
            nsegm=nsegm+1
            iresinc=0
            nbfacsumtot=nbfacsumtot+nbfacsum
            bfacsum=0.d0
            nbfacsum=0
          end if
          if (lineinp(1:6) .eq. 'HETATM' .and. iqha .eq. 0
     -        .and. ipredict .eq. 0) then
            call askyn('Do you want to include heteroatoms',34,1,1,iha,
     -        16,0)
            iqha=1
          end if
          if (lineinp(1:4) .eq. 'ATOM' .or.
     -             lineinp(1:6) .eq. 'HETATM') then
            if (lineinp(iccol1:iccol2) .eq.
     -        '************************') then
              if (iundef .eq. 0) then
                write (6,2032)
                call askyn('Do you want to drop the atom',28,1,-1,
     -            idrop,00,0)
                if (idrop .eq. 1) then
                  iundef=1
                else
                  iundef=2
                end if
              end if
              if (iundef .eq. 1) then
                lineinp(idcol:idcol)=asterisk
                nudel=nudel+1
c               print *,'n,nudel,idcol=',n,nudel,idcol
                go to 210
              end if
            end if
210         n=n+1
            if (n .eq. 1) segidprev='****'
            call checkdim(n,MAXREC,'MAXREC',6,'number of atoms',15,0)
            index(n)=nlines
            icont=0
            if (n .eq. 1 .and. nconfig .eq. 0) then
c             Check if PDB type given corresponds to chainid input
              s4=lineinp(iisegcol(1,iocpdb):iisegcol(2,iocpdb))
              s1=lineinp(iisegcol(1,iobpdb):iisegcol(2,iobpdb))
              if (s1 .eq. ' ' .and. s4 .eq. '    ') then
                if (iruntyp .ne. iruntypclean .and.
     -              nconfig .le. maxrepconf) then
                  print *,'WARNING: segment ID is missing - you may ',
     -              'have to run a CLEAN operation'
                  nosegid=1
                end if
              else if (inpcrdtyp .eq. iobpdb) then
                if (s1 .eq. ' ') then
                  write (6,2115) 'Charmm',' ','SSSS','S','    '
                  icont=1
                end if
              else
                if (s4 .eq. '    ') then
                  write (6,2115) 'Brookhaven','S','    ',' ','SSSS'
                  icont=1
                end if
              end if
            end if
            if (icont .gt. 0) call askstop(1)
            if (lineinp(1:6) .eq. 'HETATM' .and. iha .eq. 0) then
              lineinp(idcol:idcol)=asterisk
              nhdel=nhdel+1
            else
              if (lineinp(7:11) .eq. '*****') then
                iatseqnum=iatseqnumprev+1
                write (lineinp(7:11),1007) mod(iatseqnum,100000)
              else
                read (lineinp(7:11),*) iatseqnum
              end if
              if (inpcrdtyp .eq. iocpdb .or. inpcrdtyp .eq. ioa3pdb .or.
     -            inpcrdtyp .eq. ioa4pdb)
     -          call leftadjustline(lineinp,irescol1,irescol2)
              if (n .gt. 1) then
                if (lineinp(isegcol1:isegcol2) .ne.
     -              segidprev(1:nsegmid)) then
                  bfacavg(nsegm)=0.0
                  nsegm=isegno(n-1)+1
                  if (nbfacsum .gt. 0) bfacavg(nsegm)=bfacsum/nbfacsum
                  nbfacsumtot=nbfacsumtot+nbfacsum
                  bfacsum=0.d0
                  nbfacsum=0
                else if (iatseqnum .lt. iatseqnumprev .and.
     -                   iatseqnumprev .ne. 99999 .and.
     -                   ignoresegnojump .eq. 0) then
                  write (6,1216) iatseqnum,iatseqnumprev
                  if (ipredict .eq. 0) then
                    call askyn('Do you want to start a new segment',34,
     -                1,1,newsega,000,0)
                    if (newsega .eq. 0) then
                      call askyn(
     -                  'Do you want to ignore all sequence # jumps',42,
     -                  1,1,ignoresegnojump,000,0)
                      print *,'NOte, that the <C>lean option can fix ',
     -                  'such sequence number jumps'
                    end if
                  end if
                  if (nsegm .eq. isegno(n-1) .and. newsega .eq. 1) then
                    nsegm=nsegm+1
                    iresinc=0
                  end if
                end if
              end if
              isegno(n)=nsegm
              iatseqnumprev=iatseqnum
            end if
            if (lineinp(iccol1:iccol2) .eq.
     -        '************************') then
              do k=1,3
                c(k,n)=999.9
              end do
              write (lineinp(iccol1:iccol2),1102) (c(k,n),k=1,3)
            else
              read (lineinp(iccol1:iccol2),1102,err=9911) (c(k,n),k=1,3)
            end if
c           Check col 67-76 for blank
            nnbl=0
            do ic=67,min0(76,lcinp)
              if (lineinp(ic:ic) .ne. ' ') nnbl=nnbl+1
            end do
            if (nnbl .gt. 0 .and. iblank67_76 .eq. 1) then
              if (iask67_76 .eq. 1) then
                print *,'Columns 67-76 are not blank'
                call askyn('Do you want to blank out columns 67-76',38,
     -            1,1,iblank67_76,0,0)
                iask67_76=0
              end if
              if (iblank67_76 .eq. 1) call blankout(lineinp,67,76)
            end if
c           Check col 77-78 for chemical symbol
            iatnum(n)=99
            if (lineinp(77:78) .eq. '  ' .and. iread77_78 .eq. 1) then
              if (iwarn77_78 .eq. 1) then
                write (6,2031) lineinp(1:min0(79,lcinp))
                iwarn77_78=0
              end if
            else
              if (iread77_78 .eq. 0 .and. iask77_78 .eq. 1) then
                idef=1
                if (lineinp(77:78) .eq. '  ') idef=-1
                call askyn(
     -            'Do you want to read chemical symbols from col 77-78',
     -            51,1,idef,iread77_78,0,0)
                iask77_78=0
                if (n .gt. 1 .and. iread77_78 .eq. 1) write (6,2033)
                if (iread77_78 .eq. 1) nocalcium=0
              end if
              if (iread77_78 .eq. 1) then
                if (lineinp(77:77) .ne. ' ') then
                  call uplow(lineinp(78:78),lineinp(78:78),1,noabc)
                  iatnum(n)=ianum(lineinp(77:78),1,2)
                else
                  atsymbol(1:1)=lineinp(78:78)
                  atsymbol(2:2)=' '
                  iatnum(n)=ianum(atsymbol,1,2)
                end if
                if (iatnum(n) .eq. 99) write (6,2034) n,lineinp(77:78),
     -            lineinp(1:min0(79,lcinp))
              end if
            end if
            if (iatnum(n) .eq. 99)
     -        iatnum(n)=ianum(lineinp(inamcol1:inamcol2),1,nnamcol)
            call readint(lineinp,iresncol1,iresncol2,iresno(n),2,0,
     -        ierr)
            if (ierr .gt. 0) go to 9919
            if (iresno(n) .lt. iresnoprev .and.
     -        lineinp(isegcol1:isegcol2) .eq. segidprev(1:nsegmid)) then
              if (newsegr .eq. -1) then
                write (6,1217) iresno(n),iresnoprev,iatseqnum,n
                idefans=1
                if (iresnoprev .eq. 9999) idefans=-1
                call askyn('Do you want to start a new segment',34,
     -            1,idefans,newsegr,000,0)
              end if
              if (newsegr .eq. 1) then
c               Increment segment number, but only when it was not already done
                isegdec=0
                if (n .gt. 1) then
                  if (isegno(n) .gt. isegno(n-1)) isegdec=1
                end if
                if (isegdec .eq. 0) then
                  nsegm=nsegm+1
                  isegno(n)=nsegm
                end if
              end if
            end if
            if (ialtcol .gt. 0) altcol(n)=lineinp(ialtcol:ialtcol)
            if (iinscol .gt. 0) then
              inscol(n)=lineinp(iinscol:iinscol)
              if (inscol(n) .ne. ' ') ninsres=ninsres+1
            end if
            if (iresnoprev-iresno(n) .eq. incresno .and. incresno .ne.0)
     -        iresinc=iresinc+incresno+1
            segidprev(1:nsegmid)=lineinp(isegcol1:isegcol2)
            froccin(n)=1.0
            iresnoprev=iresno(n)
            iresno(n)=iresno(n)+iresinc
            if (iocccol2 .ge. iocccol1)
     -        read (lineinp(iocccol1:iocccol2),*,err=9912,end=8805)
     -        froccin(n)
            bfac=0.0
            read (lineinp(61:66),*,end=8805,err=8805) bfac
            bfacsum=bfacsum+bfac
            nbfacsum=nbfacsum+1
8805        if (iqcol2 .ge. iqcol1)
     -       read (lineinp(iqcol1:iqcol2),*,err=9912,end=8803) charge(n)
            go to 8804
c           Blank charge field should be zero
8803        charge(n)=0.0
          end if
8804      if (lineinp(1:4) .eq. 'ATOM' .or.
     -        lineinp(1:6) .eq. 'HETATM') then
            call altcolcheck(lineinp,idcol,iruntyp,ialtcol,
     -        ikeepfullalt,altnam,naltnam,naltrec,naltdel,ipredict,
     -        froccin,outfile,altfile,namleno,namlena,ncol,n,asterisk,
     -        maxrec)
            line(lineread)=lineinp
          end if
          lineprev=lineinp
        end do
23      if (nhdel .gt. 0) print *,nhdel,' heteroatoms will be deleted'
        if (nudel .gt. 0)
     -    print *,nudel,' undetermined atoms will be deleted'
        nrecdel=nhdel+nudel
        call altdelcheck(n,naltnam,naltrec,altnam,ipredict,iruntyp,
     -    altcol,ialtcol,idcol,naltdel,ncol,line,index,idrop,linewr,
     -    altfile,namlena,asterisk,maxrec)
c       Get last segment bfac average
        bfacavg(nsegm)=0.0
        if (nbfacsum .gt. 0) bfacavg(nsegm)=bfacsum/nbfacsum
        nbfacsumtot=nbfacsumtot+nbfacsum
c       Check if residue number position is ok (Charmm sometimes uses col27)
        ndig27=0
        n9999=0
        ictest=27
        do ia=1,n
          if (idigit(line(index(ia))(ictest:ictest),1) .eq. 1) then
c           Residue number extends to col 27 - fix it
            call readint(line(index(ia)),iresncol1,iresncol2+1,
     -        iresno(ia),2,0,ierr)
            if (ierr .gt. 0) go to 9911
            if (iresno(ia) .le. 9999) then
              line(index(ia))(iresncol1:iresncol2)=
     -          line(index(ia))(iresncol1+1:iresncol2+1)
              line(index(ia))(iresncol2+1:iresncol2+1)=' '
            else
              n9999=n9999+1
            end if
            ndig27=ndig27+1
c            write (77,6544) ia,iresno(ia),iresncol1,iresncol2,irno_old
c6544        format(' ia,iresno=',2i6,' ic1,2=',2i3,' irno_old=',i5)
          end if
        end do
        if (ndig27 .gt. 0) write (6,2116) ictest,ndig27,ictest
        if (n9999 .gt. 0) write (6,2117) ictest
        if (naltdel .gt. 0) then
          print *,naltdel,' alternate atoms will be deleted'
          if (ipredict .eq. 1) print *,'Predictable run keeps the ',
     -      'records marked with the 1st altternative character'
        end if
        nrecdel=nrecdel+naltdel
c        write (77,9671) (ia,line(index(ia))(iseqncol1:iseqncol2),
c     -    isegno(ia),iresno(ia),ia=1,n)
c9671    format(i6,1x,a,' isegno=',i3,' iresno=',i5)
      else if (inpcrdtyp .eq. iommod) then
c-------Macromodel structure file input
c       Decode number of atoms
        lineread=1
        read (inpt,1001,end=9999) lineinp
        line(1)=lineinp
        call lastchar(lineinp,lc,ncol)
        if (lc .le. 1 .and. lineinp(1:1) .eq. ' ') go to 9917
        read (lineinp(2:6),1007,err=9915,end=9999) n
        title=lineinp(8:87)
        trtitle(1)=lineinp(8:87)
        ntitltr=1
c       Make sure title is not blank
        icol=7
        call nextchar(lineinp,icol,132)
        if (icol .gt. 80) line(1)(7:16)='Macromodel'
        if (nconfig .le. maxrepconf) then
          call lastchar(lineinp,ifc,ncol)
          write (6,1010) lineinp(7:ifc)
        end if
        ntitlin=1
        nlines=n+1
        do i=1,n
          read (inpt,1001,end=3000) lineinp
          line(i+1)=lineinp
          do j=1,nnoblank3
            k=noblank3(j)
            if (lineinp(k:k) .ne. ' ') write (6,2051) i,k,lineinp(k:k)
          end do
          index(i)=i+1
          lineread=i+1
          read (lineinp(iccol1:iccol2),1103,err=9911) (c(k,i),k=1,3)
          read(lineinp(iqcol1:iqcol2),*) charge(i)
          call readint(lineinp,ipotcol1,ipotcol2,ityp,4,1,irerr)
          if (ityp .lt. 1 .or. ityp .gt. 64) then
            print *,'ERROR: invalid Macromodel atom type:',ityp
            iatnum(i)=99
          else
            iatnum(i)=mmatno(ityp)
          end if
        end do
c       Collect segment id-s (they can be scattered)
        nsegm=0
        do ia=1,n
          isegno(ia)=0
          s1=line(ia+1)(isegcol1:isegcol2)
          do j=1,nsegm
            if (s1 .eq. segid4(j)(1:1)) isegno(ia)=j
          end do
          if (isegno(ia) .eq. 0) then
            nsegm=nsegm+1
            segid4(nsegm)(1:1)=s1
            isegno(ia)=nsegm
          end if
          call readint(line(ia+1),iresncol1,iresncol2,iresno(ia),2,1,
     -      irerr)
        end do
      else if (inpcrdtyp .eq. iommc .or. inpcrdtyp .eq. iommc4) then
c-------MMC .slt file input
        ntitlin=0
        nsegm=1
        n=0
        do irec=1,maxrec
          lstch=0
          do while (lstch .le. 1)
            read (inpt,1001,end=41) lineinp
            call lastchar(lineinp,lstch,80)
          end do
          line(irec)=lineinp
          lineread=irec
          icol=1
          call nextchar(lineinp,icol,132)
          if (lineinp(icol:icol) .eq. '!') then
c           Comment line - see if it is title
            if (irec .eq. ntitlin+1) then
              ntitlin=ntitlin+1
              if (nconfig .eq. 0 .and. ntitlin .eq. 1)
     -          title=lineinp(1:80)
              if (nconfig .le. maxrepconf) then
                call lastchar(lineinp,ifc,ncol)
                write (6,1010) lineinp(1:ifc)
              end if
            end if
          else
            n=n+1
            index(n)=irec
            call readint(lineinp,iresncol1,iresncol2,iresno(n),2,1,
     -        irerr)
            if (n .gt. 1) then
              if (iresno(n) .lt. iresno(n-1)) nsegm=nsegm+1
            end if
            isegno(n)=nsegm
            read (lineinp(iccol1:iccol2),1101,err=9911) (c(k,n),k=1,3)
            iatnum(n)=ianum(lineinp(inamcol1:inamcol2),1,nnamcol)
            read(lineinp(iqcol1:iqcol2),*) charge(n)
          end if
        end do
        itrunc=maxrec
41      if (iclone .eq. -1) call askyn(
     -    'Do you have cloning information',31,1,-1,iclone,17,0)
        if (iclone .gt. 0) then
          write (6,*) 'NOTE: cloned atoms will have the coordinates ',
     -      'of the atom it was cloned from'
          neednnlist=1
          call getint('Number of clones',16,0,1,1000,nclone,98)
          do ic=1,nclone
42          write (question(6:8),1003) ic
            call getintline(question,51,1,n,in12,2,00)
            iaclnf(ic)=in12(1)
            iaclnl(ic)=in12(2)
            if (iaclnf(ic) .gt. iaclnl(ic))then
              print *,'Invalid range'
              go to 42
            end if
            if (ic .gt. 1) then
              if (iaclnf(ic) .lt. iaclnl(ic-1)) then
                print *,'Cloned molecules should be specified in ',
     -            'increasing order'
                go to 41
              end if
            end if
            call getint('Number of copies',16,0,1,0,ncopcln(ic),0)
          end do
          do ic=nclone,1,-1
            incr0=iaclnl(ic)-iaclnf(ic)+1
            incr=(ncopcln(ic)-1)*incr0
            incri0=index(iaclnl(ic))-index(iaclnf(ic))+1
            incri=(ncopcln(ic)-1)*incri0
c           Make room
            do irec=index(n),index(iaclnl(ic))+1,-1
              line(irec+incri)=line(irec)
            end do
c           Adjust index, segment and  and residue number
            incrres=iresno(iaclnl(ic))-iresno(iaclnf(ic))+1
            incrseg=isegno(iaclnl(ic))-isegno(iaclnf(ic))+1
            if (iaclnf(ic) .gt. 1) then
              if (isegno(iaclnf(ic)) .eq. isegno(iaclnf(ic)-1)) then
c               If cloned segment was not recognized as such, increment isegno
                do ia=iaclnf(ic),iaclnl(ic)
                  isegno(ia)=isegno(ia)+1
                end do
              end if
            end if
            do ia=n,iaclnl(ic)+1,-1
              index(ia+incr)=index(ia)+incri
              isegno(ia+incr)=isegno(ia-incr0)+ncopcln(ic)*incrseg
              iresno(ia+incr)=iresno(ia)+(ncopcln(ic)-1)*incrres
              write (line(index(ia+incr))(iresncol1:iresncol2),1007)
     -          iresno(ia+incr)
            end do
c           write (77,4411) 'B',(index(i),i=1,n+incr)
c           Clone
            do id=1,ncopcln(ic)-1
              do ir=index(iaclnf(ic)),index(iaclnl(ic))
                line(ir+id*incri0)=line(ir)
              end do
              do ia=iaclnf(ic),iaclnl(ic)
                isegno(ia+id*incr0)=isegno(ia)+id
                index(ia+id*incr0)=index(ia)+id*incri0
                iresno(ia+id*incr0)=iresno(ia)+id*incrres
                write (line(index(ia+id*incr0))(iresncol1:iresncol2),
     -            1007) iresno(ia+id*incr0)
              end do
            end do
            n=n+incr
            lineread=lineread+incri
          end do
c         write (77,4411) 'C',(index(i),i=1,n)
          do ia=1,n
            irec=index(ia)
            read (line(irec)(iccol1:iccol2),1101,err=9911)
     -        (c(k,ia),k=1,3)
c           write (77,*) 'ia,irec,an,iat=',ia,irec,
c    -       line(irec)(inamcol1:inamcol2),iatnum(ia)
c           iatnum(ia)=ianum(line(irec)(inamcol1:inamcol2),1,nnamcol)
          end do
        end if
        nlines=lineread
        if (nlines .eq. 0) iendfound=1
c        do ia=1,n
c          write (77,7622) ia,index(ia),line(index(ia))(1:80),
c     -      iresno(ia),isegno(ia)
c7622  format(2i5,a,' ires,isegno=',2i5)
c        end do
c        do ia=1,n
c          write (77,7633) ia,index(ia),iatnum(ia),(c(k,ia),k=1,3)
c7633  format(2i5,' iatnum=',i5,' c=',3f10.5)
c        end do
      else if (inpcrdtyp .eq. iogro) then
c-------Gromos/Gromacs .gro file input
        ntitlin=1
        read (inpt,1001,end=9999) lineinp
        line(1)=lineinp
        title=lineinp(1:80)
        trtitle(1)=lineinp(1:80)
        ntitltr=1
        if (nconfig .le. maxrepconf) then
          call lastchar(lineinp,ifc,ncol)
          write (6,1010) lineinp(1:ifc)
        end if
        read (inpt,1001) lineinp
        line(2)=lineinp
        call lastchar(lineinp,lc,ncol)
        if (lc .le. 1 .and. lineinp(1:1) .eq. ' ') go to 9917
        read (lineinp(1:5),1007,err=9915) n
        nlines=2
        call checkdim(n+nlines,MAXREC,'MAXREC',6,'number of lines',15,0)
        do i=1,n
          nlines=nlines+1
          index(i)=nlines
          lineread=nlines
          read (inpt,1001,end=3000) lineinp
          line(nlines)=lineinp
          read (lineinp(iccol1:iccol2),1102,err=9911) (c(k,i),k=1,3)
          call readint(lineinp,iresncol1,iresncol2,iresno(i),2,1,irerr)
          if (iqcol2 .ge. iqcol1)
     -      read(lineinp(iqcol1:iqcol2),*) charge(i)
        end do
      else if (inpcrdtyp .eq. iomol2) then
c-------Tripos .mol2 file input
        if (nconfig .eq. 0)
     -    print *,'NOTE: This format is mostly for input only'
        call read_mol2(inpt,line,nlines,n,iresno,iatnum,title,ititle,
     -    nneig,ineig,iseqncol1,iseqncol2,inamcol1,irescol1,iresncol1,
     -    iresncol2,iccol1,iccol2,ipotcol1,iqcol1,iqcol2,ncol,c,charge,
     -    index,iout,nerr,maxng,MAXREC)
        if (nerr .gt. 0) go to 9919
        innlistread=1
      else if (inpcrdtyp .eq. iomae) then
c-------Schrodinger Maestro input
        if (nconfig .eq. 0)
     -    print *,'NOTE: This format is for input only'
        call read_mae_mol(n,nbonds,iatnum,iresno,isegno,nsegm,c,charge,
     -    nneig,ineig,index,line,inpcrdtyp,iofull,inpt,
     -    iout,ierr,MAXREC,MAXNEIG)
        if (ierr .gt. 0) go to 9919
        nlines=n
        innlistread=1
      else if (inpcrdtyp .eq. iocif) then
c-------PDBx/mmCIF input
        if (nconfig .eq. 0)
     -    print *,'NOTE: This format is for input only'
        call read_cif_mol(n,iatnum,iresno,isegno,nsegm,c,froccin,charge,
     -    altcol,inscol,index,line,linefread,title,ltitle,inpcrdtyp,
     -    iofull,idcol,iruntyp,ialtcol,iinscol,ikeepfullalt,altnam,
     -    naltrec,naltdel,ninsres,ipredict,outfile,altfile,namleno,
     -    namlena,pdbid,ncol,asterisk,inpt,ierr,MAXREC)
        nrecdel=nrecdel+naltdel
        if (ierr .gt. 0) stop
        lineread=n
        innlistread=1
      else if (inpcrdtyp .eq. ioxxx) then
c-------NOT USED  input
        print *,'Input type not implemented'
        stop
      else if (inpcrdtyp .eq. ioins) then
c-------Insight .car file input
        ntitlin=4
        read (inpt,1001) lineinp
        line(1)=lineinp
c       Check for older version
        icarvers=3
        if (lineinp(17:17) .eq. '2') then
          icarvers=2
          idcol=65
          itcol=62
          iqcol1=69
          iqcol2=75
          iresncol1=57
          iresncol2=61
        end if
        print *,'Version',icarvers
        read (inpt,1001) lineinp
        line(2)=lineinp
        call writeline(6,lineinp,1,ncol,0)
        ipbc=0
        if (lineinp(1:6) .eq. 'PBC=ON' .or. lineinp(1:6) .eq. 'PBC=2D')
     -    ipbc=1
        read (inpt,1001) lineinp
        line(3)=lineinp
        title=lineinp(1:80)
        trtitle(1)=lineinp(1:80)
        ntitltr=1
        if (nconfig .le. maxrepconf) then
          call lastchar(lineinp,ifc,ncol)
          write (6,1010) lineinp(1:ifc)
        end if
        read (inpt,1001) lineinp
        line(5)=lineinp
        if (ipbc .eq. 1) then
          ntitlin=ntitlin+1
          read (inpt,1001) lineinp
          line(ntitlin)=lineinp
        end if
        n=0
        lineread=ntitlin
        ntitadd=0
        do i=1,maxrec
          lstch=0
          do while (lstch .le. 1)
            read (inpt,1001,end=51) lineinp
            call lastchar(lineinp,lstch,80)
          end do
          line(lineread+1)=lineinp
          lineread=lineread+1
          call checkdim(lineread,MAXREC,'MAXREC',6,'number of lines',15,
     -      0)
          if (lineinp(1:1) .eq. '!') then
            if (lineread .eq. ntitlin+i) then
c             Additional header lines
              ntitadd=ntitadd+1
              call writeline(6,line(lineread),1,ncol,0)
            end if
          else if (lineinp(1:3) .eq. 'end') then
            nsegm=nsegm+1
          else
            n=n+1
            index(n)=lineread
            if (n .gt. 1) then
              if (lineinp(isegcol1:isegcol2) .ne.
     -            lineprev(isegcol1:isegcol2)) nsegm=nsegm+1
            end if
            isegno(n)=nsegm
            lineread=ntitlin+i
            read (lineinp(iccol1:iccol2),1104,err=9911) (c(k,n),k=1,3)
            atomnam='     '
            atomnam(1:2)=lineinp(ichemcol1:ichemcol2)
            iatnum(n)=ianum(atomnam,1,2)
            icol=iresncol1
            call nextchar(lineinp,icol,132)
            call nextblank(lineinp,icol,132)
            icoll=icol-1
            do while (icoll .gt. iresncol1 .and.
     -                idigit(lineinp(icoll:icoll),1) .eq. 0)
              icoll=icoll-1
            end do
            call readint(lineinp,iresncol1,icoll,iresno(n),2,1,irerr)
            read(lineinp(iqcol1:iqcol2),*) charge(n)
          end if
          lineprev=lineinp
        end do
51      ntitlin=ntitlin+ntitadd
        if (icarvers .eq. 3) then
          rnu='      '
          nrescar=0
          do i=1,n
            if (line(index(i))(iresncol1:iresncol2) .ne. rnu) then
              nrescar=nrescar+1
              rnu=line(index(i))(iresncol1:iresncol2)
              write (line(index(i))(iresncol1:iresncol2),1006) nrescar
            else
              write (line(index(i))(iresncol1:iresncol2),1006) nrescar
            end if
          end do
          write (6,2028) nrescar
        end if
      else if (inpcrdtyp .eq. ionxyz .or. inpcrdtyp .eq. iosxyz .or.
     -         inpcrdtyp .eq. iosxyzrq) then
c-------Insight free format input
c       Decode number of atoms
        call lastchar(lineinp,lc,ncol)
        if (lc .le. 1 .and. lineinp(1:1) .eq. ' ') go to 9917
        read (inpt,*,err=9915) n
        ntitlin=0
        call checkdim(n,MAXREC,'MAXREC',6,'number of atoms',15,0)
        do i=1,n
          read (inpt,1001,end=3000) lineinp
          line(ntitlin+1)=lineinp
          index(i)=ntitlin+i
          lineread=ntitlin+i
          if (inpcrdtyp .eq. ionxyz) then
            read (lineinp,*,err=9913) iatnum(i),(c(k,i),k=1,3)
          else if (inpcrdtyp .eq. iosxyz .or.
     -             inpcrdtyp .eq. iosxyzrq) then
            icol=1
            call nextchar(lineinp,icol,132)
            i1=icol
            call nextblank(lineinp,icol,132)
            atsymbol='  '
            atsymbol(1:icol-i1)=lineinp(i1:icol-1)
            atomnam='     '
            atomnam(1:2)=atsymbol
            iatnum(i)=ianum(atomnam,1,2)
            if (inpcrdtyp .eq. iosxyz) then
              read (lineinp(icol:ncol-icol+1),*,err=9911) (c(k,i),k=1,3)
            else
              read (lineinp(icol:ncol-icol+1),*,err=9916)
     -          (c(k,i),k=1,3),iresno(i),charge(i)
            end if
          end if
        end do
      else
c-------Grasp .crg file
        ntitlin=1
        read (inpt,1001) lineinp
        line(ntitlin)=lineinp
        do while (line(ntitlin)(1:1) .eq. '!')
          ntitlin=ntitlin+1
          read (inpt,1001) lineinp
          line(ntitlin)=lineinp
        end do
        ntitlin=ntitlin-1
        n=1
        lineread=ntitlin
        do while (.true.)
          read (inpt,1001,end=3001) lineinp
          line(ntitlin+n)=lineinp
          index(n)=ntitlin+n
          lineread=ntitlin+n
          n=n+1
        end do
3001    n=n-1
      end if
      go to 9900
3000  itrunc=i-1
      go to 9919
9999  iendfound=1
      if (ioktoend .eq. 1) return
      go to 9919
9911  write (6,9001) lineinp(iccol1:iccol2),iccol1,iccol2,icinc
      print *,'n=',n
      go to 9919
9912  print *,'Invalid syntax for charge'
      go to 9919
9913  print *,'Invalid syntax for coordinates or atomno'
      go to 9919
9915  print *,'Invalid syntax for number of atoms'
      go to 9919
9916  print *,'Invalid syntax for coordinates or resno or charge'
      go to 9919
9917  print *,'Number of atoms is missing'
c     stop
9919  if (lineread .eq. 0) then
        write (6,1210)
      else
        write (6,1209) lineread,lineinp
      end if
      if (iendfound .eq. 1 .or. n .eq. 0) then
        if (nconfig .eq. 0) then
          print *,'ERROR: Coordinate file ',inpfile(1:namleni),
     -    ' is empty'
          iaskagain=1
        else
          print *,'WARNING: run out of data after ',nconfig,
     -      ' configurations'
          iendfound=1
          if (nconfig .eq. 0) call askstop(0)
        end if
      end if
      if (itrunc .gt. 0) then
        print *,'WARNING: Coordinate file contains only ',itrunc,
     -    ' atoms instead of ',n
        n=itrunc
        if (nconfig .eq. 0) call askstop(0)
        itrunc=0
      end if
      if (iunrecog .eq. 1) then
        print *,'Contents of file ',inpfile(1:namleni),' are not in ',
     -    formatname(inpcrdtyp),' format'
        if (nconfig .eq. 0) then
          iaskagain=1
        else
          stop
        end if
      end if
      if (nconfig .gt. 0 .and. iendfound .eq. 0) then
        print *,'Error occurred at molecule number ',nconfig
c       Skip to the next molecule for some of the input formats
        iskipok=0
        if (ispdb(inpcrdtyp) .eq. 1) then
c         PDB
          do while (lineinp(1:3) .ne. 'END')
            call blankout(lineinp,1,ncol)
            read (inpt,1010,end=9900) lineinp
          end do
          iskipok=1
        else if (ischarmm(inpcrdtyp) .eq. 1) then
c         Charmm
          do while (lineinp(1:1) .ne. '*')
            call blankout(lineinp,1,ncol)
            read (inpt,1010,end=9900) lineinp
          end do
          linefread=1
          iskipok=1
        else if (inpcrdtyp .eq. iomol2) then
c         Tripos mol2
          do while (lineinp(1:17) .ne. '@<TRIPOS>MOLECULE')
            call blankout(lineinp,1,ncol)
            read (inpt,1010,end=9900) lineinp
          end do
          linefread=1
          iskipok=1
        else if (inpcrdtyp .eq. iomae) then
c         Schroinger Maestro
          do while (lineinp(1:6) .ne. 'f_m_ct')
            call blankout(lineinp,1,ncol)
            read (inpt,1010,end=9900) lineinp
          end do
          linefread=1
          iskipok=1
        end if
        if (iskipok .eq. 1) then
          print *,'Skipped input to the next molecule'
          go to 100
        end if
      end if
9900  if (naltnam .gt. 0) close (40)
      if (nbfacsumtot .eq. 0) call zeroit(bfacavg,nsegm)
      return
1001  format(a132)
1003  format(i3)
1006  format(i4)
1007  format(i5)
1010  format(a)
1101  format(3f10.5)
1102  format(3f8.3)
1103  format(3f12.5)
1104  format(3f15.9)
1105  format(3f20.10)
1015  format(1x,a)
1213  format(' ERROR: number of title lines (',i9,') exceeds 32')
1215  format(' ERROR: File ',a,' does not',/, 8x,'appear to be in ',a,
     -  ' format')
1216  format(' Sequence number',i6,' is less than the previous',
     -  ' sequence number (',i6,')')
1217  format(' Residue id ',i6,' is less than the previous residue id ',
     -  '(',i6,')',/,' Sequence # read:',i7,' n=',i7,')')
1209  format(' Last line read (',i7,'-th):',/,a)
1210  format(' No data found in the file')
2028  format(' Residue sequence names are changed to sequence numbers',
     -  /,' Number of residues found=',i4)
2031  format(' WARNING: column 77-78 is blank - atomic number will be ',
     -  'deduced from atom name:',/,1x,a)
2032  format('Lines with undefined coordinates were found',/,
     -  ' The atom may be dropped or the coordinates may be set to ',
     -  '999.9')
2033  format(' WARNING: Not all ATOM/HETATM records have chemical ',
     -  'symbols')
2034  format(' ERROR: invalid chemical symbol found in record # ',i6,
     -  ':',a,/,1x,a)
2051  format(' WARNING: non-blank character in line of atom ',i5,
     -  ' column ',i3,':',a1)
2077  format(' Note: all heteroatoms will be kept and',/,
     -  7x,'only the first of alternate records will be used')
2115  format(' WARNING: Input structure has ',a,'-type segment ID, ',
     -  'i.e.,'/,'ATOM  99999 AAAA RRR ',a1,
     -  ' 999       0.0     0.0     0.0    1.00  0.0       ',a4,/,
     -  ' instead of',/,'ATOM  99999 AAAA RRR ',a1,
     -  ' 999       0.0     0.0     0.0    1.00  0.0       ',a4)
2116  format(' Column ',i2,' is not blank for',i7,' atoms',/,
     -  5x,'- residue numbers will be read including column',i3,/,
     -  5x,'- residue numbers < 10000 will be shifted 1 colum to the ',
     -  'right')
2117  format(' NOTE: residue numbers > 9999 will be left using ',
     -  'column',i3)
9001  format(' Invalid syntax for coordinates:',a,
     -  ' (cols ',i3,' - ',i3,'; icinc=',i1,')')
      end
