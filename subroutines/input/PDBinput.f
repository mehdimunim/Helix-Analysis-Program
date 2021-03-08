PDB input
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
c