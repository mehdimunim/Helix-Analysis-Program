module simulaid
contains
          character*(*) inpfile,outfile
          character* 132 line(maxrec)
          character*80 title,trtitle(32),molname,linewr
          character*4 segid4(maxrsd),pdbid
          character*1 altcol(maxrec),inscol(maxrec)
          common /clonedat/ nclone,iaclnf(1000),iaclnl(1000),ncopcln(1000)
          character*11 formatname
          common /formats/ iqdconv(20),formatname(19)
          character*5 crdext
          common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
          character*2 iatnm2
          common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
          common /askCA/ iaskca,nocalcium
          common /logging/ logfile,ipredict
          character*1 altnam(50),s1,asterisk,ans
          character*2 atsymbol
          character*4 s4,segidprev
          character*8 atomnam
          character*8 segmid
          character*6 rnu
          character*51 question
          character* 132 lineinp,lineprev,linenext
          character*200 altfile
          call setcol(inpcrdtyp,ncol,idcol,ialtcol,iinscol,
    c-------Charmm input
            call blankout(lineinp,1,ncol)
                call checkforetot(1,lineinp,nconfig,etotread,ietotread,1)
                call lastchar(lineinp,ifc,ncol)
    c       Decode number of atoms
            call lastchar(lineinp,lc,ncol)
            call checkdim(nlines,MAXREC,'MAXREC',6,'number of lines',15,0)
              call readint(lineinp,iresncol1+icinc,iresncol2+icinc,
    c         Blank charge field should be zero
                call leftadjustline(lineinp,isegcol1,isegcol2)
    c-------PDB input
                call askyn(
              call blankout(lineinp,1,ncol)
              call lastchar(lineinp,lcinp,ncol)
              call checkforetot(6,lineinp,nconfig,etotread,ietotread,1)
                call nextblank(lineinp,ic,ncol)
                  call blankout(lineinp,1,ncol)
                call checkdim(nlines,MAXREC,'MAXREC',6,'number of lines',15,
    c           Don't check for data after END for trajectory combining
    c             See if file ends
                  call blankout(linenext,1,ncol)
                  call lastchar(linenext,ilc,ncol)
                  call askyn(
    c             Replace END with TER
                call askyn('Do you want to include heteroatoms',34,1,1,iha,
                    call askyn('Do you want to drop the atom',28,1,-1,
    c               print *,'n,nudel,idcol=',n,nudel,idcol
                call checkdim(n,MAXREC,'MAXREC',6,'number of atoms',15,0)
    c             Check if PDB type given corresponds to chainid input
                        call askyn('Do you want to start a new segment',34,
                          call askyn(
                    c(k,n)=999.9
    c           Check col 67-76 for blank
                    call askyn('Do you want to blank out columns 67-76',38,
    c           Check col 77-78 for chemical symbol
                    call askyn(
                      call uplow(lineinp(78:78),lineinp(78:78),1,noabc)
                call readint(lineinp,iresncol1,iresncol2,iresno(n),2,0,
                    call askyn('Do you want to start a new segment',34,
    c               Increment segment number, but only when it was not already done
    c           Blank charge field should be zero
                call altcolcheck(lineinp,idcol,iruntyp,ialtcol,
            call altdelcheck(n,naltnam,naltrec,altnam,ipredict,iruntyp,
    c       Get last segment bfac average
    c       Check if residue number position is ok (Charmm sometimes uses col27)
    c           Residue number extends to col 27 - fix it
                call readint(line(index(ia)),iresncol1,iresncol2+1,
    c            write (77,6544) ia,iresno(ia),iresncol1,iresncol2,irno_old
    c6544        format(' ia,iresno=',2i6,' ic1,2=',2i3,' irno_old=',i5)
    c        write (77,9671) (ia,line(index(ia))(iseqncol1:iseqncol2),
    c     -    isegno(ia),iresno(ia),ia=1,n)
    c9671    format(i6,1x,a,' isegno=',i3,' iresno=',i5)
    c-------Macromodel structure file input
    c       Decode number of atoms
            call lastchar(lineinp,lc,ncol)
    c       Make sure title is not blank
            call nextchar(lineinp,icol,132)
              call lastchar(lineinp,ifc,ncol)
              call readint(lineinp,ipotcol1,ipotcol2,ityp,4,1,irerr)
    c       Collect segment id-s (they can be scattered)
              call readint(line(ia+1),iresncol1,iresncol2,iresno(ia),2,1,
    c-------MMC .slt file input
                call lastchar(lineinp,lstch,80)
              call nextchar(lineinp,icol,132)
    c           Comment line - see if it is title
                    call lastchar(lineinp,ifc,ncol)
                call readint(lineinp,iresncol1,iresncol2,iresno(n),2,1,
              call getint('Number of clones',16,0,1,1000,nclone,98)
                call getintline(question,51,1,n,in12,2,00)
                call getint('Number of copies',16,0,1,0,ncopcln(ic),0)
    c           Make room
    c           Adjust index, segment and  and residue number
    c               If cloned segment was not recognized as such, increment isegno
    c           write (77,4411) 'B',(index(i),i=1,n+incr)
    c           Clone
    c         write (77,4411) 'C',(index(i),i=1,n)
    c           write (77,*) 'ia,irec,an,iat=',ia,irec,
    c    -       line(irec)(inamcol1:inamcol2),iatnum(ia)
    c           iatnum(ia)=ianum(line(irec)(inamcol1:inamcol2),1,nnamcol)
    c        do ia=1,n
    c          write (77,7622) ia,index(ia),line(index(ia))(1:80),
    c     -      iresno(ia),isegno(ia)
    c7622  format(2i5,a,' ires,isegno=',2i5)
    c        end do
    c        do ia=1,n
    c          write (77,7633) ia,index(ia),iatnum(ia),(c(k,ia),k=1,3)
    c7633  format(2i5,' iatnum=',i5,' c=',3f10.5)
    c        end do
    c-------Gromos/Gromacs .gro file input
              call lastchar(lineinp,ifc,ncol)
            call lastchar(lineinp,lc,ncol)
            call checkdim(n+nlines,MAXREC,'MAXREC',6,'number of lines',15,0)
              call readint(lineinp,iresncol1,iresncol2,iresno(i),2,1,irerr)
    c-------Tripos .mol2 file input
            call read_mol2(inpt,line,nlines,n,iresno,iatnum,title,ititle,
    c-------Schrodinger Maestro input
            call read_mae_mol(n,nbonds,iatnum,iresno,isegno,nsegm,c,charge,
    c-------PDBx/mmCIF input
            call read_cif_mol(n,iatnum,iresno,isegno,nsegm,c,froccin,charge,
    c-------NOT USED  input
    c-------Insight .car file input
    c       Check for older version
            call writeline(6,lineinp,1,ncol,0)
              call lastchar(lineinp,ifc,ncol)
                call lastchar(lineinp,lstch,80)
              call checkdim(lineread,MAXREC,'MAXREC',6,'number of lines',15,
    c             Additional header lines
                  call writeline(6,line(lineread),1,ncol,0)
                call nextchar(lineinp,icol,132)
                call nextblank(lineinp,icol,132)
                call readint(lineinp,iresncol1,icoll,iresno(n),2,1,irerr)
    c-------Insight free format input
    c       Decode number of atoms
            call lastchar(lineinp,lc,ncol)
            call checkdim(n,MAXREC,'MAXREC',6,'number of atoms',15,0)
                call nextchar(lineinp,icol,132)
                call nextblank(lineinp,icol,132)
    c-------Grasp .crg file
    c     stop
    c       Skip to the next molecule for some of the input formats
    c         PDB
                call blankout(lineinp,1,ncol)
    c         Charmm
                call blankout(lineinp,1,ncol)
    c         Tripos mol2
                call blankout(lineinp,1,ncol)
    c         Schroinger Maestro
                call blankout(lineinp,1,ncol)
          character*1 altnam(maxrec),asterisk
          character*132 lineinp
          character*200 outfile,altfile
    c     print *,'ALTCOLCHECK n=',n,' ialtcol=',ialtcol
    c         Alternative location found
                  call askyn(
                  call changeext(outfile,altfile,namleno,namlena,
                  call openfile(40,0,' ',1,'new',altfile,namlena,
    c           See if name already occurred
                call writeline(40,lineinp,1,ncol,0)
    c           Duplicates will be dropped
          character*1 altnam(50),altcol(maxrec),asterisk
          character*80 linewr
          character*132 line(maxrec)
          character*200 altfile
    c     print *,'ALTDELCHECK n=',n,' ialtcol=',ialtcol,
    c    -  ' ipredict=',ipredict
              call askyn(linewr,39,1,0,idrop,0,0)
              call askstop(0)
            call askyn(
    c         Mark locations to be deleted
          character*1 ans
          character*8 atomnam
          character*80 line
          character*200 filename
          call quiz(ans,ians,'n',' ',0,'charge input',12,0,5,6,0)
            call askyn('Do you have charges for the solvents too',40,1,1,
            call openfile(inpt,0,'Amber prtmtop',13,'old',filename,namlen,
              charge(i)=charge(i)/18.2223
    c         write (77,*) i,' q=',charge(i)
            call openfile(inpt,0,'Charmm PSF',10,'old',filename,namlen,
            call find_n_psf(inpt,line,nerr,nread,natspsf,'!NATOM',6)
              call blankout(line,1,80)
                call nextstring(line,ic,ic1,ic2,80)
              charge(ia)=0.0
            call openfile(inpt,0,'Autodock .pdbqt',15,'old',filename,namlen,
                charge(ia)=0.0
            call askstop(1)
          character*80 line
          character*(*) label
            call blankout(line,1,80)
            call nextstring(line,ic,ic1,ic2,80)
            call nextstring(line,ic,ic21,ic22,80)
          character*80 linewr
    c     print *,'WRITEBOND iout_bond,iout_conn=',iout_bond,iout_conn,
    c    -  ' nslt=',nslt,' maxng=',maxng
    c         write (6,8693) ia,(ineig(j,ia),j=1,nneig(ia))
    c8693     format(i4,' nn:',10i4)
          character*132 line(maxat)
    c       print *,line(il)(1:77)
    c       print *,'nlfound=',nlfound
          call shiftmol(c,n,shift,c,1.0)
          call rotate_c(c,n,rot,c,'REVERSE',7)
          character*(*) inpfile,outfile
          character*132 line(maxat),blankline
          character*4 segid4(maxseg)
          character*8 atnames(maxat),resnames(maxrsd)
          character*2 iatnm2
          common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
          character*1 abc,digits,hexdigits
          common /charactersets/ ihex(25),abc(62),digits(14),hexdigits(25)
          character*1 xyz
          common /axislab/ xyz(3)
          character*1 segids(62),chainprev,chlist(20),chidused(60),chnew
          character*2 nm2
          character*132 linew
          character*200 cellfile,lineinp
    c      write (6,9782) 'A',(segid4(i),i=1,nsegslt)
    c9782  format(1x,a,' SEGID4:',20(a4,'|'))
          call zeroit(abcabg,6)
    c       Make atom names conform PDB convention
          call askyn(
            cellfile=inpfile
            call openfile(inptcell,0,'cell information',16,'old',cellfile,
    c     Find crytal info first
                call nextstring(lineinp,ic,ic1,ic2,200)
                call nextstring(lineinp,ic,ic1,ic2,200)
                call nextstring(lineinp,ic,ic1,ic2,200)
                call nextstring(lineinp,ic,ic1,ic2,200)
                call nextstring(lineinp,ic,ic1,ic2,200)
                call nextstring(lineinp,ic,ic1,ic2,200)
          call unitmat(uxyz)
            call trnsfr(edge,abcabg,3)
            call getxyz('Edge length in the ',19,' direction (A)',14,
            call getreal('ey-ez angle(in deg)',19,90.0,abcabg(4),1,0)
            call getreal('ex-ez angle(in deg)',19,90.0,abcabg(5),1,0)
            call getreal('ex-ey angle(in deg)',19,90.0,abcabg(6),1,0)
    c       Generate cell's unit vectors
            cosxy=cos(abcabg(6)/radtodeg)
              cellxyz(k,i)=edge(i)*uxyz(k,i)
            call checkdim(n*nasym,maxat,'MAXREC',6,
    c       Now create the transformed coordinates
              call rotate_c(c,n,rotasym(1,1,isym),c(1,(isym-1)*n+1),
              call shiftmol(c(1,(isym-1)*n+1),n,shiftasym(1,isym),
    c       Generate the unit cell vertices
            call zeroit(vertex,3*8)
            call compact_ucell(c,ctemp,itemp,molsltlim,ntot,nasym*nsegslt,
    c       icompact=0
    c       if (irect .eq. 1) then
    c         call try_to_compact(c,ctemp,molsltlim,edge,n,nasym*nsegslt,
    c    -      icompact)
    c         if (icompact .eq. 1) then
    c           print *,'Further translations made the cell more compact:'
    c           call extension(c,itemp,0,1,n,cmin,cmax,c0,1,0,v)
    c         end if
    c       else
    c         print *,'NOTE: the Edit/trans... option may translate the ',
    c    -      'images to a more compact form'
    c       end if
              call extension(c,itemp,0,1,ntot,cmin,cmax,c0,1,0,vol)
    c         Gather all contacts
                call cofms(c(1,(is-1)*n+1),crm_as(1,is),n,atw)
              call trnsfr(ctemp,c,n*3)
                      call zeroit(shift,3)
    c                 Test copy # is shifted by shift
                      call arrsum(crm_as(1,is),shift,crmtest,3)
    c                 if (dist2(crm_as(1,1),crm_test) .lt. d_org*1.2) then
    c                   See if this copy is in contact with asym #1
                          call arrsum(c(1,(is-1)*n+ia),shift,ctest,3)
    c                         Use this copy
                              call checkdim(n*ncopy,maxat,'MAXREC',6,
                              call trnsfr(ctemp(1,(ncopy-1)*n+1),
                              call shiftmol(ctemp(1,(ncopy-1)*n+1),n,shift,
    c                 end if
              call trnsfr(c,ctemp,3*n*ncopy)
              call lastchar(line(i),lc,80)
              chainprev=' '
                  chainprev=linew(22:22)
                call lastchar(linew,lc,80)
                call askyn(
                  call zeroit(shift,3)
                    call getint(
                  call shiftmol(c,nw,shift,ctemp,1.0)
                    chainprev=' '
                        chainprev=linew(22:22)
                      call lastchar(linew,lc,80)
              call askyn('Do you want to write the cell vertices/edges too',
    c       Generate biological oligomers
              chidused(ic)=segid4(ic)(1:1)
            call zeroiti(iusedseg,0,nsegslt)
                  call blankout(lineinp,1,200)
                  call lastchar(lineinp,lc,200)
                  call nextchar(lineinp,ic,200)
                    chlist(nchread)=lineinp(ic:ic)
                    call blankout(lineinp,1,200)
    c         Read CIF info
                call blankout(lineinp,1,200)
              call lastchar(lineinp,lc,200)
              call nextchar(lineinp,ic,200)
                chlist(nchread)=lineinp(ic:ic)
    c         print *,'CHLIST=',(chlist(i),i=1,nchread)
                call blankout(lineinp,1,200)
              call lastchar(lineinp,lc,200)
                call blankout(lineinp,1,200)
                call lastchar(lineinp,lc,200)
                  call nextstring(lineinp,ic,ic1,ic2,200)
    c             print *,'I=',i,' IC1,2=',ic1,ic2
                    call blankout(lineinp,1,200)
                    call lastchar(lineinp,lc,200)
                call blankout(lineinp,1,200)
    c         Check if transformation causes actual change
    c           No chain was specified - use all chains
    c           Set mask for chains selected
                call zeroiti(maskseg,0,nsegslt)
    c         print *,'MASKSEG=',(maskseg(i),i=1,nsegslt)
    c               Generate new chain ID
                    chnew=abc(ich)
                    chidused(nsegtot+1)=chnew
                    chnew=segid4(ic)(1:1)
    c              write (6,9671) (chidused(i),i=1,nsegtot)
    c9671          format(' CHIDUSED=',20a1)
    c             Transform chain ic and write the corresponding PDB file
                  call rotate_c(c(1,is1o),natss,rotasym(1,1,isym),
                  call shiftmol(ctemp(1,isn),natss,shiftasym(1,isym),
    c             print *,'IC =',ic ,' if,l=',(molsltlim(k,ic),k=1,2)
    c             print *,'ICN=',icn,' if,l=',(molsltlim(k,icn),k=1,2)
    c               Just change the coordinates and chain ids
                      call lastchar(linew,lc,80)
    c               Create full new record
                      call blankout(linew,1,132)
                      call createrec_fromcif(linew,iocif,iotyp,
                      call lastchar(linew,lc,132)
          character*132 linew,blankline
          character*1 segnam1
          character*4 atnam,segnam,chemnam
          character*6 potnam
          character*8 resnam
          character*2 iatnm2
          common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
    c     print *,'CREATEREC_FROMCIF iat=',iat,' IRES=',ires,' IAN=',iatnum
          chemnam(1:2)=iatnm2(iatnum)
          call createrec(linew,iocif,iotyp,ctemp(1),ctemp(2),ctemp(3),' ',
    c     subroutine createrec(line,inpcrdtyp,ioutyp,cx,cy,cz,
    c    -  atnam,resnam,segnam,iseqno,iresnum,iresid,chemnam,potnam,
    c    -  frocc,q,nqdec,iqspace,mmcgm,nn,in,nconfig,ibnd,ihetat,blankline)
    c     character* 132 line,blankline
    c     character*1 resnam1
    c     character*4 atnam,segnam,chemnam
    c     character*6 potnam
    c     character*8 resnam
    c     character*2 mmcgm
    c         print *,'IFST,ILST=',ifst,ilst
              call trnsfr(ctemp,c(1,ifst),3*nats)
              call trnsfr(ctemp(1,nats+1),vertex,24)
              call extension(ctemp,itemp,0,1,nats+8,cmin,cmax,c0,1,0,vol)
              call checkonedir(c,ctemp,itemp,ifst,ilst,edge,uxyz,
              call checkonedir(c,ctemp,itemp,ifst,ilst,edge,uxyz,
              call checkonedir(c,ctemp,itemp,ifst,ilst,edge,uxyz,
    c     print *,'COMP NSHIFT=',nshift
    c       Translate back and forth along the x axis
    c         Generate translation along the iaxis-th axis
    c         print *,'SHIFT=',shift
              call shiftmol(c(1,ifst),nats,shift,ctemp,float(ix))
              call trnsfr(ctemp(1,nats+1),vertex,24)
              call extension(ctemp,itemp,0,1,nats+8,cmin,cmax,c0,1,0,volnew)
                call trnsfr(c(1,ifst),ctemp,3*nats)
              call countout_rect(ifst,ilst,c,edge,nouttot,kmax,nslt)
    c           Translate
                call zeroit(shift,3)
                call shiftmol(c(1,ifst),ilst-ifst+1,shift,ctemp(1,ifst),
    c           Check if ctemp is more compact
                call countout_rect(ifst,ilst,ctemp,edge,nouttotnew,kmax,
                  call trnsfr(c(1,ifst),ctemp(1,ifst),3*(ilst-ifst+1))
          call zeroiti(nxyzp,0,3)
          call zeroiti(nxyzm,0,3)
    c      write (6,1000) nxyzm,nxyzp
    c      write (6,1001) ifst,ilst,nouttot,edge
    c1000  format(' COUNTOUT nxyzm=',3i9,' nxyzp=',3i9)
    c1001  format(' COUNTOUT ifst,ilst=',2i6,' nout=',i6,' e=',3f10.5)
          character*(*) title
          character* 132 line(maxat)
          character*1 xyz
          common /axislab/ xyz(3)
    c     character*2 bondord(maxng,maxat)
    c     character*4 resnamprev
          character*80 lineinp
    c     print *,'READMOL2 iomol2,maxng,maxat=',iomol2,maxng,maxat
    c     print *,'READMOL2 inamcol1=',inamcol1,' irescol1=',irescol1,
    c    -  ' iresncol1=',iresncol1
    c     resnamprev='XYZU'
          call blankout(lineinp,1,80)
            call blankout(lineinp,1,80)
    c     print *,line(nlines_mol+2)(1:60)
          call lastchar(line(nlines_mol+1),ltitle,80)
    c     print *,'nats,nbonds=',nats,nbonds
            call blankout(lineinp,1,80)
              call blankout(line(nlines),1,80)
              call blankout(line(nlines),1,80)
    c         Atom id
              call nextchar(lineinp,ic,lenrec)
              call nextblank(lineinp,ic,lenrec)
    c         Atom name
              call nextchar(lineinp,ic,lenrec)
              call nextblank(lineinp,ic,lenrec)
    c         Check for blank within the atom name
              call nextchar(lineinp,icc,lenrec)
    c           Fill in the blank
    c         Coordinates
                call nextchar(lineinp,ic,lenrec)
                call nextblank(lineinp,ic,lenrec)
    c         Atomtype
              call nextchar(lineinp,ic,lenrec)
              call nextblank(lineinp,ic,lenrec)
    c         Substance id (residue number)
              call nextchar(lineinp,ic,lenrec)
              call nextblank(lineinp,ic,lenrec)
    c         Residue (substructure) name
              call nextchar(lineinp,ic,lenrec)
              call nextblank(lineinp,ic,lenrec)
    c         if (resnamprev(1:ic-ic1) .ne. lineinp(ic1:ic-1)) then
    c           New residue
    c           nres=nres+1
    c           resnamprev(1:ic-ic1)=lineinp(ic1:ic-1)
    c         end if
    c         Charge
              call nextchar(lineinp,ic,lenrec)
              call nextblank(lineinp,ic,lenrec)
    c         Status bit - not read for now
              call askyn('There were 25 errors so far. Do you want to stop',
            call askstop(1)
          call zeroiti(nnmol2,0,natsmol2)
            call blankout(lineinp,1,80)
    c       Bond order label is ignored for now
    c       ic=1
    c       do ifnd=1,4
    c         call nextchar(lineinp,ic,lenrec)
    c         ic1=ic
    c         call nextblank(lineinp,ic,lenrec)
    c       end do
    c         bondord(nnmol2(ia),ia)(1:ic-ic1)=lineinp(ic1:ic-1)
    c         bondord(nnmol2(ja),ja)(1:ic-ic1)=lineinp(ic1:ic-1)
          character*132 line(maxat)
          character*400 linein
          character*20 items(MAXCOL)
          character*50 colid(MAXCOL)
          character*2 iatnm2
          common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
          character*4 segn4,segn4prev
    c     print *,'READ_MAE inpt,iout,inptyp_mae=',inpt,iout,inptyp_mae
          call setcol(inptyp_mae,ncol,idcol,ialtcol,iinscol,
          call read_head_mae('m_atom',6,nats,linein,len,inpt,iout)
    c     print *,'NATS=',nats
    c     Read item descriptors
          call read_colid_mae(ndcol,colid,lcolid,'atoms',5,linein,
    c     print *,'ix_x,y,z=',ix_x,ix_y,ix_z
    c     print *,'ix_rnam,rnum,anam,atno=',ix_rnam,ix_rnum,ix_anam,ix_atno
              c(k,ia)=999.9
            call blankout(line(ia),inamcol1,inamcol2)
            call blankout(line(ia),1,ncol)
            call blankout(linein,1,len)
              call nextstring(linein,ic,ic1,ic2,len)
    c         write (6,*) 'IA,ICOL,IC1,IC2=',ia,icol,ic1,ic2
    c        write (6,8791) ia,(ic,litems(ic),items(ic)(1:litems(ic)),ic=1,3)
    c8791  format(' IA=',i2,3(' ic=',i2,' l=',i2,' item=',a))
              call blankout(line(ia),irescol1,irescol2)
    c    -      resnam(ia)(1:litems(icol))=items(icol)(1:litems(icol))
    c         if (ian(ia) .le. 99)
    c    -      line(ia)(ichemcol1:ichemcol1+1)=iatnm2(ian(ia))
    c         print *,'IA=',ia,' INCHEMCOL1=',ichemcol1,' cnam=',
    c    -      iatnm2(ian(ia))
    c         Just added hydrogen
                  call writeint(line(ia),ic,nha,nhlen)
    c         All (UNK) names are blank - generate them from atomic numbers
              call zeroiti(iancount,0,99)
                    call writeint(line(ia),ic,iancount(ian(ia)),nal)
    c     print *,'NATS=',nats
          call read_head_mae('m_bond',6,nbonds,linein,len,inpt,iout)
    c     print *,'NBONDS,maxng=',nbonds,maxng
    c       Read bond info
            call read_colid_mae(ndcol,colid,lcolid,'bonds',5,linein,
            call zeroiti(nn,0,nats)
            call zeroiti(in,0,maxng*nats)
    c       print *,'NDCOL=',ndcol
              call blankout(linein,1,len)
              call nextstring(linein,ic,ic1,ic2,len)
              call nextstring(linein,ic,ic1,ic2,len)
              call nextstring(linein,ic,ic1,ic2,len)
    c     print *,'icol=',icol,' litems=',litems(icol)
          character*(*) lab,line
            call blankout(line,1,len)
            call nextchar(line,ic,len)
    c           Header found
                call findnextchar(']',line,ic,len)
          character*(*) lab,line
          character*50 colid(maxcol)
    c     print *,'READ COLID inpt,iout=',inpt,iout,' lab=',lab(1:llab)
          colid(1)='index'
            call blankout(line,1,200)
            call lastchar(line,lc,len)
              call nextchar(line,ic,len)
    c           print *,'ndcol,ic,lc=',ndcol,ic,lc
                colid(ndcol)(1:lcolid(ndcol))=line(ic:lc)
          character*1 altcol(maxat),inscol(maxat),altnam(50),asterisk
          character*132 line(maxat),linewr
          character*80 title
          character*200 linein,outfile,altfile
          character*4 pdbid
          character*2 iatnm2,an2
          common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
          character*4 chain,chain_prev
          character*8 atomname
          character*10 rectyp
    c     print *,'READ_CIF_MOL inptyp_cif=',inptyp_cif
          call setcol(inptyp_cif,ncol,idcol,ialtcol,iinscol,
          call get_title_cif(inpt,linein,lineread,title,ltitle,pdbid)
          call askyn('Do you want to use the author (PDB) residue numbers',
    c     print *,'TITLE:',title(1:ltitle)
          call zeroiti(ic_xyz,0,3)
            call nextblank(linein,ic,200)
            call blankout(linein,1,200)
          chain_prev='    '
            call blankout(line(nats),1,132)
              c(k,nats)=999.9
            call lastchar(linein,lc,200)
              call nextstring(linein,ic,limcol(1,nc),limcol(2,nc),200)
              call save_cif_rec(rectyp,10,'iiccol',6,limcol(1,icr),
                call save_cif_rec(rectyp,10,'iiresncol',9,
    c           Residue number has to be generated - temporarily set to zero
                call save_cif_rec(rectyp,10,'iiresncol',9,
              call save_cif_rec('Atom name',9,'iinamcol',8,
              call save_cif_rec('Residue name',12,'iirescol',8,
              call save_cif_rec('Segment id',10,'iisegcol',8,
              chain='    '
              chain(1:limcol(2,icr)-limcol(1,icr)+1)=
                chain_prev=chain
                call uplow(an2(2:2),an2(2:2),1,noabc)
              call findname(an2,iatnm2,1,99,ian(nats),2)
              call save_cif_rec('Chemical sym',12,'iichemcol',9,
              call save_cif_rec(rectyp,9,'ialtcol',7,limcol(1,icr),
              call save_cif_rec(rectyp,9,'iinscol',7,limcol(1,icr),
    c     subroutine save_cif_rec(label,llabel,clabel,lclabel,limcol,
    c    -  line,linein,ic1,ic2,len,noadjust,ierr)
              call save_cif_rec(rectyp,9,'iiocccol',8,limcol(1,icr),
              call save_cif_rec(rectyp,4,'iiqcol',6,limcol(1,icr),
              call save_cif_rec(rectyp,10,'iiseqncol',9,limcol(1,icr),
    c         Generate residue number
    c           New residue
            call altcolcheck(line(nats),idcol,iruntyp,ialtcol,
    c       write (77,*) nats,' C=',(c(k,nats),k=1,3)
    c       write (77,*) nats,' IRESNO=',iresno(nats)
    c       write (77,*) nats,' ATNAME=',line(nats)(inamcol1:inamcol2)
    c       write (77,*) nats,' RESNAME=',line(nats)(irescol1:irescol2)
    c       write (77,*) nats,' IATNO=',ian(nats)
    c       write (77,*) nats,' ISEGNUM=',isegno(nats)
    c       write (77,*) nats,' ISEQNUM=',ii
    c       write (77,*) nats,' OCCPANCY=',frocc(nats)
    c       write (77,*) nats,' BETA=',charge(nats)
    c       write (77,*) nats,' BETA(L)=',line(nats)(iqcol1:iqcol2)
    c       write (77,*) nats,' SEGID=',line(nats)(isegcol1:isegcol2)
    c       write (78,*) nats,' ',line(nats)
    c       print *,nats,'----------------------------------'
            call blankout(linein,1,200)
            call askstop(0)
          call altdelcheck(nats,naltnam,naltrec,altnam,ipredict,iruntyp,
    c     do i=1,nats
    c       write (78,7861) nats,naltrec,ialtcol,altcol(i),line(i)
    c     end do
    c7861 format(i6,' naltrec=',i4,' ialtcol=',i3,' altcol=',a,/,
    c    -  ' line=',a,/)
          character*80 title,titles(4)
          character*200 linein
          character*4 pdbid
    c     rewind inpt
          call zeroiti(len,0,4)
              call blankout(linein,1,200)
              call blankout(linein,1,200)
              call lastchar(linein,lc,200)
              call nextstring(linein,ic,ic1,ic2,200)
              call nextstring(linein,ic,ic1,ic2,200)
          call getint('Choice (1-4):',12,2,1,4,in,000)
          call blankout(title,1,80)
            call getname(title,ltitle,'New title',9,80,'TITLE',5,1,000,00)
          character*(*) label,clabel
          character*132 line
          character*200 linein
    c     print *,'SAVE_CIF_REC ',label(1:llabel),' ic1,i2=',ic1,ic2
    c     print *,' REC:',linein(limcol(1):limcol(2)),' LIM:',limcol
    c     if (ic1 .le. 111 .and. ic2 .ge. 111 .and.
    c    -  label(1:llabel) .ne. 'Alt at id ') then
    c       print *,label(1:llabel),' ic1,ic2=',ic1,ic2
    c       stop
    c     end if
          character* 132 line(maxrec)
          character*8 resnames(maxrsd)
          character*8 ionresnam(100),molresnam(100)
    c     Check for ions - make them separate molecules
          call askyn('Do you have ions in this system',31,1,-1,ions,18,0)
          call askyn('Do you have molecular residues in this system',
    c       Gather residue names
    c       Last block of solute residues may be ions
    c        write (77,7272) (is,(molsltlim(k,is),k=1,2),is=1,nmolslt)
    c7272    format(i4,' molsltlims=',2i6)
    c       Make requested residues molecules
    c             Do nothing
    c
    c             Last residue of the segment
    c             First residue of the segment
    c             Residue in the middle of the segment
    c       write (77,7272) (is,(molsltlim(k,is),k=1,2),is=1,nmolslt)
          character*4 anames(n)
          character*4 typnames(200),number
            call findname(anames(ia),typnames,1,ntyps,ix,4)
          call zeroiti(ntyp,0,200)
            call findname(anames(ia),typnames,1,ntyps,ix,4)
              call lastchar(anames(ia),lc,4)
              call writeint(number,ic,ntyp(ix),len)
          character*4 segid4
          character*8 segid8
          character* 132 line
    c     Get a max 4 character segment name
            call leftadjustline(segid8,1,nsegcol)
    c     Calculate molar volume using known values of amino and nucleic acids
          character* 132 line(maxrec)
          character*8 resname,resnamela
          character*1 aanames1
          character*2 mmodtoamb
          character*3 aanames3
          common /atnamcon/ mmodtoamb(100),aanames1(58),aanames3(58),
          common /pmvol/ pmvaana(58)
    c     print *,'MOLARV ir1,ir2,irn1,irn2=',ir1,ir2,irn1,irn2
    c     print *,'MOLARV nslt=',nslt
            call readint(line(index(ia)),irn1,irn2,iresnum,2,1,irerr)
    c         New residue starts
              call leftadjustn(resname,resnamela,lresn)
              call findname(resnamela,aanames3,1,naanames,ix,3)
                call findname(resnamela,aanames3,
    c         write (77,*) 'ia=',ia,' resname=',resname(1:3),
    c    -      ' naa,nnw,nna,nnf=',naa,nnw,nna,nnf
            call leftadjustn(resname,resnamela,lresn)
          character*1 dssplab(maxrsd)
          character*1 charl1(50),charl2(50)
          character*8 atnam
          character* 132 line(maxrec)
          character*1 typc
          character*21 ssname
          common /dsspnames/ lssname(9),ssname(9),typc(9)
    c     print *,'DSSP n1,n,nslt,iwdssp=', n1,n,nslt,iwdssp
              call trnsfr(ch(1,nres+1),c(1,ia),3)
    c         New residue
    c           print *,'IRESO,N=',ireso,iresn
    c           print *,'ICFOUND,INFOUND,IOFOUND,IHFOUND=',
    c    -          icfound,infound,iofound,ihfound
    c             Generate H coordinates from C(prev), CA and N
                  call zeroit(ch(1,nres),3)
                        ch(k,nres)=c(k,ixn(nres))+rx(k)/sqrt(dnh)
    c                 print *,'nconfig,maxrepconf=',nconfig,maxrepconf
                  call trnsfr(ccprev,c(1,ixc(nres)),3)
                  call zeroit(ccprev,3)
    c       write (77,*) 'ia,ireso,nres=',ia,ireso,nres
    c      do i=1,nres
    c        write (77,1144) i,'CA:',line(index(ixa(i)))(1:80)
    c        write (77,1144) i,'C :',line(index(ixc(i)))(1:80)
    c        write (77,1144) i,'N :',line(index(ixn(i)))(1:80)
    c        write (77,1144) i,'O :',line(index(ixo(i)))(1:80)
    c        dd=dist2(ch(1,i),c(1,ixn(i)))
    c        ddco=dist2(c(1,ixc(i)),c(1,ixo(i)))
    c        ddcn=dist2(c(1,ixc(i)),c(1,ixn(i)))
    c        write (77,*) 'dCO=',sqrt(ddco),' dCN=',sqrt(ddcn)
    c        write (77,1133) i,sqrt(dd),(ch(k,i),k=1,3)
    c1133    format(i5,' d(N-H)=',f10.5,' ch=',3f10.5)
    c1144    format(i6,1x,a,1x,a)
    c1145    format(i6,1x,a,1x,3f10.5)
    c      end do
              call trnsfr(cres(1,ir),c(1,ixc(ir)),3)
              call zeroit(cres(1,ir),3)
          call nnlistsim(1,nres,cres,nneig,ineig,indices,nbox,
    c                if (eij .lt. threshold .or. eji .lt. threshold)
    c     -            write (77,1156) ir,jr,eij*(0.42*0.20*332.0),
    c     -              eji*(0.42*0.20*332.0)
    c1156            format(' ir,jr=',2i5,' eij,eji=',2f10.5)
    c                  if (eij .lt. threshold .or. eji .lt. threshold)
    c                  if (ir .lt. 999)
    c     -              write (77,1155) ir,jr,eij,eji,threshold
    c1155              format(' ir,jr=',2i5,' eij,eji=',2f10.5,' tr=',f10.5)
    c                  if (ir .lt. 999) write (77,1154) ir,ihbneig(ir),
    c     -              enghb(ir),jr,ihbneig(jr),enghb(jr)
    c1154              format(' ir,ihbneig(ir),enghb(ir)=',2i5,f10.5,
    c     -              ' jr,ihbneig(jr),enghb(jr)=',2i5,f10.5)
    c     A hydrogen bond exist from C=O(i) to NH(ihbneig(i))
          call zeroiti(iparal,0,nres)
          call zeroiti(iantiparal,0,nres)
    c       Look for the next SS element
    c       write (77,*) 'Start check ir=',ir,' nhbinc=',nhbinc
    c         Possible helix start
    c           write (77,6622) ir,nhbinc1,nhbinc2,nhbinc3,ihfound
    c6622       format(' ir=',i4,' dihn1,2,3=',3i3,' ifhound=',i2)
    c           write (06,*) 'Helix nss,ifss,ilss,itypss=',nss,ifss(nss),
    c    -        ilss(nss),itypss(nss)
    c         Possible Lambda helix
    c           write (77,6622) ir,nhbinc1,nhbinc2,nhbinc3,ihfound
    c6622       format(' ir=',i4,' dihn1,2,3=',3i3,' ifhound=',i2)
    c           write (77,*) 'L Helix nss,ifss,ilss,itypss=',nss,ifss(nss),
    c         Possible sheet
    c           write (77,*) 'Sheet? ir,jr,ihbneig(jr)=',ir,jr,ihbneig(jr)
    c               HB(i,j) = HB(j,i) => Antiparallel
    c               write (77,*) 'A1 ir=',ir
    c                 HB(j-1,i)=HB(i,j+1) => Parallel
    c                 write (77,*) 'P2 ir=',ir
    c                 HB(i-1,j)=HB(j,i+1) => Parallel
    c                 write (77,*) 'P1 ir=',ir
    c                 HB(i-1,j+1)=HB(j-1,i+1) => Antiparallel
    c                 write (77,*) 'A2 ir=',ir
    c           write (77,*)'End do isfound,ir=',isfound,ir
    c         write (77,*)'islast,isfirst,npar,napar,ineigmin,max=',
    c    -       islast,isfirst,npar,napar,ineigmin,ineigmax
    c           write (77,*) 'Sheet nss,ifss,ilss,itypss=',nss,ifss(nss),
    c    -        ilss(nss),itypss(nss)
    c         Neither helix nor sheet - just skip over
    c         do while (ir .lt. nres .and. ihbneig(ir) .ne. 0)
    c         end do
    c        write (77,1001) iss,ifss(iss),ilss(iss),itypss(iss)
    c1001    format(' SS',i4,' Start at',i4,' End at',i4,' Type=',i2)
                charl1(ic)=' '
                charl2(ic)=' '
                  call angles(dist2(c(1,ixa(ir-2)),c(1,ixa(ir))),
                  cosa=dble(cbend)
    c             write (77,*) ir,' bend=',bend
                  charl1(ir-iresf+1)='?'
    c             write (77,*) ir,' tors=',tors*radtodeg
                    charl2(ir-iresf+1)='+'
                    charl2(ir-iresf+1)='-'
                  charl2(ir-iresf+1)='?'
    c       Tentative alternative for turn detection (see Rose's 1977 paper)
            call normplane(c(1,ixa(1)),c(1,ixa(3)),c(1,ixa(5)),rnprev)
            call normplane(c(1,ixa(2)),c(1,ixa(3)),c(1,ixa(4)),rn1prev)
              call radcirc(c(1,ixa(ir-2)),c(1,ixa(ir)),c(1,ixa(ir+2)),r)
              call normplane(c(1,ixa(ir-2)),c(1,ixa(ir)),c(1,ixa(ir+2)),rn)
              call normplane(c(1,ixa(ir-1)),c(1,ixa(ir)),c(1,ixa(ir+1)),rn1)
              call angles(dist2(c(1,ixa(ir-2)),c(1,ixa(ir))),
              cosa=dble(cbend)
              charl1(1)=' '
              call trnsfr(rnprev,rn,3)
              call trnsfr(rn1prev,rn1,3)
          character*1 dssplab(maxrsd)
          character*6 hxoklab(3)
          character*(*) lab
          common /logging/ logfile,ipredict
          character*200 listfile
          character*80 listprompt,line
    c     Get a list of integers (nlist=1) or a pair of integers (nlist=2)
    c     from a file or from the terminal
            call askyn(
          call getname(listfile,llistfile,'Name of the list file',21,80,
            call openfile(listinp,1,'list',4,'old',listfile,llistfile,
            call readfreelist(list,len,listinp,1,maxlen)
    c       Input list interactively
              call blankout(line,1,80)
              call getname(line,lline,'0/r/#:',6,80,'',0,1,0,0)
    c           Input a range
                call getrange(ifst,999999,ilst,999999,increment,1,'number',
    c           if (irangeinp .eq. 1) then
    c             len=len+1
    c             list(1,len)=ics
    c             listprev=ics
    c           end if
    c       Input pair list
              call blankout(line,1,80)
                call getname(line,lline,listprompt,28,80,'',0,1,0,0)
    c       Create pair list from input list
              call blankout(line,1,80)
                call getname(line,lline,listprompt,31,80,'',0,1,0,0)
          character*80 line
              call blankout(line,1,80)
              call lastchar(line,lc,80)
            call nextstring(line,ic,ic1,ic2,80)
          character*39 pairprompt
          character*80 line
          character*200 outfiletmp
    c     print *,'GETCLUSTERPAIRS idebug,maxpairs,maxclustermem=',
    c    -  idebug,maxpairs,maxclustermem
          call openfile(listinp,1,'cluster member list',19,'old',outfiletmp,
              call getint(pairprompt,39,0,1,maxclustermem,np1,0)
              call getname(line,lline,'Members (comma-separated):',26,80,
              call getint(pairprompt,39,1,1,maxclustermem,np2,0)
              call getname(line,lline,'Members (comma-separated):',26,80,
          character*(*) list(maxlen),label
          character*80 line
          character*200 outfiletmp
          common /logging/ logfile,ipredict
            call blankout(line,1,80)
            call lastchar(line,iclast,80)
    c         Switch to reading from a file
              call openfile(listinp,1,'name list',9,'old',outfiletmp,
                call blankout(line,1,80)
                call lastchar(line,iclast,80)
                call blankout(list(nlist),1,lench)
                call askstop(0)
                call blankout(list(nlist),1,lench)
          character*11 formatname
          common /formats/ iqdconv(20),formatname(19)
          common /columnlim/ incol(19),iidcol(19),iialtcol(19),iiinscol(19),
          character*1 ans
            call quiz(ans,ians,'d',' ',0,'file format',11,0,5,6,0)
          character*(*) filex
          character*5 crdext
          character*(*) inoutlab
          common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
          character*11 formatname
          common /formats/ iqdconv(20),formatname(19)
          character*1 ans
          character*5 ext
          character*132 line
    c     print *,'GETCRDTYP filex=',filex(1:lenfile)
    c       Extract extension
    c       Specify PDB type
    c       Check file for segid/chaind
    c         File exists, find out PDB version
                call quiz(ans,ians,'b',' ',0,'PDB file type',13,0,5,6,0)
    c         New file - Brookhaven is assumed
    c         Determine Charmm CRD type
              call blankout(line,1,132)
              call lastchar(line,lc,132)
              call askyn('Is that OK',10,1,1,ians,00,0)
              call quiz(ans,ians,'o',' ',0,'Charmm CRD file type',20,0,5,6,
    c       Specify MMC version
            call askyn('Is the .slt file in the old format',34,1,-1,is4,134,
    c         Macromodel/Xcluster
          character*(*) analfile
          character*4 extinp
          character*6 extinp6
    c     Remove coordinate extension
          character*8 resnamslv
          character*5 crdext
          common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
          character*8 resnamsv
          common /names/ resnamsv(19)
          character*1 ansfrm
          character*6 inoutlab
            call quiz(ansfrm,ians,' ',inoutlab,6,'trajectory file format',
            call quiz(ansfrm,ians,' ',' ',0,'MMC trajectory file format',26,
          character* 132 line(maxrec),blankline
          character*80 title
          character*1 altcol(maxrec),inscol(maxrec)
          character*4 segnames(maxrsd),pflsv(100)
          character*8 atnames(maxrec),resnames(maxrsd),namesv(100)
          character*6 marker(16)
          common /columnlim/ incol(19),iidcol(19),iialtcol(19),iiinscol(19),
          character*11 formatname
          common /formats/ iqdconv(20),formatname(19)
          character*5 crdext
          common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
          character*2 iatnm2
          common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),mmatno(64),
          character*1 abc,digits,hexdigits
          common /charactersets/ ihex(25),abc(62),digits(14),hexdigits(25)
          common /logging/ logfile,ipredict
          character*2 mmcgm
          character*4 segnam,chnam,xseg
          character*6 potnam
          character*8 resnam,resnamslv,atomnam
          character*8 segnamlong
          character*31 question
    c     print *,'WRITECONF nslt=',nslt,' n0,n=',n0,n,
    c    -  ' iwriteatsym=',iwriteatsym
          call set_pdbinfo(iotyp,iwriteatsym,iwritecon,iobpdb,iocpdb,0)
          call setcol(inpcrdtyp,ncol,idcol,ialtcol,iinscol,
    c       Macromodel
            call nnlist(nslt,islvw,naslv,n,iatnum,ifchrg,c,nneig,nneiga,
            call bondord(iatnum,mmtype,n,nneig,ineig,nhneig,ibnd,maxng,c,
    c       Get segment id name list
            call getseg4(segnames(1),line(index(1)),isegcol1,nsegcol)
    c       segnames(1)(1:nsegcol)=line(index(1))(isegcol1:isegcol2)
    c         write(77,*) 'nconfig=',nconfig,' ia=',ia,'isegno=',isegno(ia)
    c         write(77,*) 'index=',index(ia)
                call getseg4(segnam,line(index(ia)),isegcol1,nsegcol)
    c           segnam(1:nsegcol)=line(index(ia))(isegcol1:isegcol2)
                call leftadjust4(segnam,segnam)
                  call lastchar(segnam,lc,4)
                    call getname(segnam,len,question,lq,4,segnam(1:lc),lc,0,
    c             For Charmm output, change the segment ID when there is room
    c             No blank, leave it
                  call leftadjustn(segnamlong,segnamlong,8)
                  call nextchar(segnamlong,ic,8)
    c         Make sure new segment ID's are different
              call zeroiti(iabc,0,62)
    c       Create segment id's from segment numbers
    c     Create atom records
          chnam='    '
          call zeroiti(ihetat,0,n)
    c         Full information is available
    c         if (iat .lt. 25) print *,'RESNAM=',resnam
    c    -      line(index(iat))(irescol1+icinc:irescol2+icinc)
              call leftadjustn(resnam,resnam,8)
              call nextchar(line(index(iat)),ic1,132)
              call nextblank(line(index(iat)),ic2,132)
    c           call readint(line(index(iat)),ic1,ic2-1,ires,2,1,irerr)
    c         Add solvent information
              charge(iat)=qsv(iatslv)
    c         Only limited information is available
            chnam(1:2)=iatnm2(iatnum(iat))
                    call askyn('Do you want to use all 6 characters',35,0,
    c          write (77,7611) iat,ires,q,segnam
    c7611      format(' iat=',i6,' ires=',i5,' q=',f10.5,' segnam=',a)
                call trnsfr(cw,c(1,iat),3)
                    cw(k)=c(k,iat)*10.0
                    cw(k)=c(k,iat)*0.1
                call createrec(line(irec),inpcrdtyp,iotyp,cw(1),cw(2),cw(3),
                call writefree(line(irec),iotyp,chnam,c(1,iat),c(2,iat),
    c       Left justify residue and atom names
              call leftadjustline(line(index(iat)),iirescol(1,iotyp),
              call rightadjustline(line(index(iat)),
              call leftadjustline(line(index(iat)),iirescol(1,iotyp),
          call writeout(iout,inpcrdtyporg,iotyp,line,index,isegno,n,marker,
    c     Puts out the generated configuration, with proper headers, separators
          character*6 marker(16)
          character*80 title,linew
          character* 132 line(maxrec),blankline,ansline,pline
          character*5 crdext
          common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
          common /logging/ logfile,ipredict
          character*11 formatname
          common /formats/ iqdconv(20),formatname(19)
          common /columnlim/ incol(19),iidcol(19),iialtcol(19),iiinscol(19),
          common /boxdat/ xtlabc(6),xtlabc0(6),box(3),box0(3),edge_gen(3,3),
          common /analparm/ nsltref_f,nsltref_l,rcut_cv,icvtyp
    c     write (6,*) 'WRITEOUT inpcrdtyp,iotyp,n,iout=',
    c    -  inpcrdtyp,iotyp,n,iout
    c         Write headers
    c           Charmm or PDB outputs
                  call nextblank(line(1),icol,132)
                  call writeline(iout,pline,1,icol+1,0)
                  call nextchar(line(1)(51:130),icol,132)
                  call lastchar(line(1),ifc,130)
                    call writeline(iout,pline,1,ifc-(icol+50)+8,0)
                    call blankout(linew,1,80)
                    call lastchar(linew,lc0,80)
                    call lastchar(title,lct,80)
    c                 Make sure not to write more than 30 lines, and no blanks
                        call blankout(linew,1,80)
                        call lastchar(linew,lc,80)
                          call lastchar(linew,lc,80)
    c                   Skip energy record since that was only in the input template
                          call nextchar(line(i),icol,132)
                          call lastchar(line(i),ifc,mxcl)
                            call writeline(iout,pline,1,ifc-icol+8,0)
                        call blankout(linew,1,80)
                        call lastchar(linew,lc,80)
                      call lastchar(line(i),ifc,mxcl)
    c                     PDB to PDB - keep full non-atom lines
                          call writeline(iout,line(i),1,ifc,0)
                          call writeline(iout,pline,1,ifc-mncl+8,0)
    c             Add data column annotation, if required
                    call writeline(iout,pline,1,0,0)
                    call writeline(iout,pline,1,0,0)
    c       Write number of atoms
                call askstop(1)
              call lastchar(ansline,ifc,incol(inpcrdtyp))
              call writeline(iout,ansline,1,ifc,0)
    c       Left-adjust residue ID
              call leftadjustline(line(index(ia)),iresidcol1,iresidcol2)
          call lastchar(line(index(1)),ifc,ncol)
            call writeline(iout,line(index(1)),1,ifc,0)
    c       REMARK records can only be placed safely if index is unchanged
              call askyn(
    c      if (keeprem .gt. 0) write (77,8811)
    c     -  (i,index(i),line(index(i))(1:80),i=1,index(n))
    c8811  format(i5,' index-',i5,' line=',a)
    c     ixprev=max0(0,index(1)-3)
    c         Segment end
                    call lastchar(line(i),lc,80)
                    call writeline(iout,line(i),1,lc,0)
            call lastchar(line(index(ia)),ifc,ncol)
            call writeline(iout,line(index(ia)),1,ifc,0)
          character*(*) line
          character*2 iatnm2
          common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
          character*2 chemsym
          character*4 name
    c     print *,'ADDATSYM ifc=',ifc
    c     Add atomic symbol to col 77-78
          chemsym=iatnm2(iatno)
             chemsym(2:2)=chemsym(1:1)
             chemsym(1:1)=' '
            call blankout(line,ifc+1,80)
          character* 132 line,blankline
          character*1 altcol,inscol,resnam1
          character*4 atnam,segnam,chemnam
          character*6 potnam
          character*8 resnam
          character*5 crdext
          character*2 mmcgm
          common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
    c      write (77,7711) ioutyp,potnam,atnam,resnam,iresnum
    c7711  format(' CREATEREC ioutyp,potnam,atnam,resnam=',i4,3(2x,a5),
    c    -   ' iresnum=',i5)
    c      write (77,*) 'CREATEREC ioutyp,iogro=',ioutyp,iogro
    c      write (80,7611) iresnum,iresid,q,segnam,altcol
    c7611  format(' CREATEREC ires,resid=',2i5,' q=',f10.5,' segnam=',a,
    c    -   ' ALTCOL=',a,'|')
    c     print *,'CREATEREC inpcrdtyp,ioutyp=',inpcrdtyp,ioutyp
          call setcol(inpcrdtyp,ncol,idcol,ialtcol,iinscol,
          call setcol(ioutyp,ncol,idcol,ialtcol,iinscol,
            call putreal(line(iqcol1+icinc:iqcol2+icinc),iqcol2-iqcol1+1,
            call leftadjustline(line,iresidcol1+icinc,iresidcol2+icinc)
            call putreal(line(iqcol1:iqcol2),iqcol2-iqcol1+1,q,2*nqdec)
            call leftadjustline(line,iresidcol1,iresidcol2)
              call putreal(line(55:60),6,frocc,0)
              call putreal(line(iqcol1+1:iqcol2),iqcol2-iqcol1,q,nqdec)
              call putreal(line(iqcol1:iqcol2),iqcol2-iqcol1+1,q,nqdec)
    c         Right shift names that are less than four characters
            call writeitem(line,ipotcol1,ipotcol2,potnam,0)
            call putreal(line(iqcol1+1:iqcol2),iqcol2-iqcol1,q,nqdec)
    c       Create also the 1-letter residue code
            call changeprot(resnam,resnam1,2)
            call writeitem(line,ipotcol1,ipotcol2,potnam,0)
            call putreal(line(iqcol1:iqcol2),iqcol2-iqcol1+1,q,nqdec)
            call writeitem(line,ichemcol1,ichemcol2,chemnam,0)
              call writeitem(line,ipotcol1,ipotcol2,potnam,0)
              call writeitem(line,ipotcol1,ipotcol2,chemnam,0)
            call putreal(line(iqcol1:iqcol2),iqcol2-iqcol1+1,q,nqdec)
          call writeitem(line,irescol1,irescol2,resnam,nrescol_i)
          call writeitem(line,isegcol1,isegcol2,segnam,min0(4,nsegcol_i))
          call writeitem(line,namcol1,namcol2,atnam,min0(4,nnamcol_i))
    c     Place one item into a record
          character* 132 line
          character*(*) name
          call blankout(line,icol1,icol2)
          common /logging/ logfile,ipredict
            call askyn('Do you want to continue',23,1,idef,icont,20,0)
          character*(*) out
          character*25 out25
    c     Puts a real number into out
    c     Check for exceeding the length
          call nextchar(out25,ic,25)
    c     print *,'PUTREAL val=',val,' out=',out,' len=',len
          character*4 chemnam
          character* 132 line
          character*5 crdext
          common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
    c     Limited free-format outputs for InsightII
          character*4 atnam
          character*3 resnam
          character*1 segnam
          character*200 analfile
    c     Write a Grasp .crg file from convdat
          character*4 ires
          character*8 convdat
          character*200 sfilename
          common /savedat/ mxresdat,maxcondat,ifst(1000),ilst(1000),
          character*200 analfile
    c     Write a UHBD .dat file from convdat
          character*4 ires
          character*8 convdat
          character*200 sfilename
          common /savedat/ mxresdat,maxcondat,ifst(1000),ilst(1000),
          character*7 libname(3)
          common /columnlim/ incol(19),iidcol(19),iialtcol(19),iiinscol(19),
    c     print *,'SETCOL inpcrdtyp=',inpcrdtyp
          character*132 inpline
    c     X-cluster new coordinates
    c     print *,'natoms,inpt,nconfig=',natoms,inpt,nconfig
    c     print *,'inpline=',inpline
          character*(*) line
          call nextchar(line,icol,132)
          character* 132 line
    c     print *,'RIGHTADJUSTLINE ic1,ic2=',ic1,ic2,' line=',line(ic1:ic2)
    c     print *,'nshift=',nshift,' line=',line(ic1:ic2)
    c*****Left-adjust a string of four characters
          character*(*) in
          character*4 out,outt
    c*****Left-adjust a string of n characters
          character*(*) in
          character*(*) out
            call blankout(out,1,n)
          common /rotmat/ matrot0(4,4),matrot(4,4),nomat0
          character*1 ans,ans0,anslc,xyzorig(3)
          common /askrotdat/ angle,ans0
          common /depthcuedat/ near,ifar,ramp0,idepth,idepthon,idrawh,
          character*1 xyz
          common /axislab/ xyz(3)
          character*3 onoff(2)
            call getreal('New angle increment',19,angle,angle,0,36)
            call uplow(ans,anslc,1,noabc)
    c     print *,'iaxis=',iaxis,' ans=',ans,' angle=',angle
          character*4 incvff,outnam,xxxx
          character*2 cvffname(100)
          character*4 charmmname(100)
    c     Convert cvff types into Charmm22 type (first approximation)
    c      write (6,8888) cvffname
    c8888  format(20(1x,a2))
          call zeroiti(iused,0,maxtyp)
          character*2 nnam
          character*4 cvffnam,iatnam,obl,ob,ct2,ct3
    c     'Refine' conversion by differentiation of some cases
    c-----Determine ring membership type
    c     print *,'REFINECONV n,maxng,maxtyp=',n,maxng,maxtyp
    c     Count ring member neighbours
    c       Set iringtyp(i) to -1 for ring atoms
    c     Set iringtyp(i) to 5 and 6, resp when ring size is known for sure
    c     Set iringtyp(i) to 4 for ring junction atoms (to be refined later)
    c     For ring nitrogens not set yet, deduce ring size from neighbours'
    c     For ring junction, set type if possible
    c     Now -1: some ring atom; 0: No ring atom; 1:5/6 junction; 2: 5/5 junction
    c     3: 6/6 junction; 4: some junction; 5:5-membered ring; 6:6-membered ring
    c         Separate carbonyl oxygen in acetic acid
          character*4 resnam,atnam(maxat),potnam(maxat)
          character*(*) outfile
    c-----Print residue name, charge
    c-----Print atoms
    c-----Print bonds
          character* 132 line(maxrec)
          character*6 marker
          character*5 crdext
          common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
    c     Rearrange coordinates to keep atoms of the same residue together
    c     print *,'SORTATRES,iresncol1,2',n,iresncol1,iresncol2
          call indexit(indexo,1,n,0)
    c     Individual intervals consist of single elements at start
          call indexit(ifa,1,n,0)
          call indexit(ila,1,n,0)
    c     Merge pairs of intervals
            call mergerec(line,index,indexo,isegno,iresncol1,iresncol2,
    c     Take care of last (odd) interval
          call mergerec(line,index,indexo,isegno,iresncol1,iresncol2,
          character*200 inpfile,numfile
          character*1 separatorchar
          common /filenuminfo/ iaskunderscore,separatorchar
    c     inout=-1: input file; inout=+1: output file; inout=+2: unpacking output f
    c     print *,'FNAMNUM inout=',inout,' iaskunderscore=',iaskunderscore,
    c    -  ' separatorchar=',separatorchar
    c     Find last '.'
    c       Skip '.1.' if present
    c         Unpacking traj or configs - always ask separator character
              call askyn('Use _ as the filenumber separator instead',
                  call askyn(
                  call askyn('Use _ as the filenumber separator instead',
    c     Insert numeral
          call writeint(numfile,nl1,n,lenw)
    c     Add part of the filename after the '.' (if any)
    c     Unpack the next conformation from the stack
          character*(*) molname
          character*132 line(maxrec)
          character*200 inpfile,outfile
          common /logging/ logfile,ipredict
    c     print *,'UNPACKCONF nextconfsel,nconfig,ionefile=',
    c    -                    nextconfsel,nconfig,ionefile
    c         Selection found - increment nextconfsel
    c         Configuration not selected - skip
    c     print *,'Unpackconf nconfig,numsel,nextconfsel=',
    c    -  nconfig,numsel,nextconfsel
    c         Currently only DOCK PDB files have molecule names read
              call filenamenum(inpfile,namleni,outfile,nl1,ifnumw,+2)
            call openfile(iuunp,0,'output',6,'new',outfile,nl1,notfnd,2,1,1,
              call askyn('Do you want to overwrite all existing files',43,
              call openfile(iuunp,0,'output',6,'new',outfile,nl1,notfnd,2,1,
    c       Open overall file (input filename with_sel added to it)
            call openfile(iuunp,0,'output',6,'new',outfile,nl1,notfnd,0,1,
            call writeline(iuunp,line(i),1,ncol,0)
          character*1 ans
          character*6 rq
          call quiz(ans,iax,' ',' ',0,'rotation type',13,0,5,6,0)
    c     print *,'ans=',ans
            call getreal(
            call unitmat(rot)
    c       print *,'ix,iy,iz=',ix,iy,iz
    c       Input rotation matrix
                call getreal(rq,6,999999.0,rot(i,j),0,00)
    c      do i=1,3
    c        write (6,7777) i,(rot(j,i),j=1,3)
    c7777    format(i3,' rot i =',3f7.4)
    c      end do
          character*4 atnam,pflab
          character* 132 line(maxrec)
          character*2 iatnm2
          common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
          character*5 crdext
          common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
          character*1 ans,modtyp
          character*2 ATD4typ
          character*80 question
          character*21 ngdefnum
          call setcol(inpcrdtyp,ncol,idcol,ialtcol,iinscol,
          chargech=0.0
            call quiz(modtyp,iansmut,' ',' ',0,'modification type',17,
              call getint('Atom number to replace',22,0,1,nslt,ia,0)
              call listatom(line,index,iatno,ia,inpcrdtyp,iofull,
              call getname(atnam,len,question,lq,4,'',0,0,0,3)
                call listatom(line,index,iatno,ia1,inpcrdtyp,iofull,
                call getint('Valence',7,nval(iatno(ia)),1,8,newn,0)
                call defaultbondl(newn,nneig(ia1),iatno(ia),iatno(ia1),
                call getreal('New distance',12,r12,rij,1,0)
                call arrdiff(c(1,ia),c(1,ia1),c12,3)
                  c(k,ia)=c(k,ia1)+rij*c12(k)/rnorm
              call getint(
              call listatom(line,index,iatno,ia1,inpcrdtyp,iofull,
              call getname(atnam,len,question,lq,4,'',0,0,0,3)
                call findroot('R2',ia1,0,iatno,nneig,ineig,line,index,
                call findroot('R3',ia2,ia1,iatno,nneig,ineig,line,index,
                  call getint('Valence',7,nval(iatnoadd),1,8,newn,0)
                  call defaultbondl(nneig(ia1)+1,newn,iatno(ia1),iatnoadd,
                  call getreal('R1-X distance',13,r12,rij,1,0)
                    call getreal('R2-R1-X angle',13,999999.0,aijk,1,0)
                    call getreal('R3-R2-R1-X torsion',18,999999.0,tijkl,0,0)
                  call addatom(1,c(1,ia3),c(1,ia3),c(1,ia2),c(1,ia1),
                call zeroiti(iadef,0,3)
    c               Neighbor exists
                    call listatom(line,index,iatno,iadef(ind),inpcrdtyp,
                    call getint('Index of the next definig atom',30,0,1,
                    call listatom(line,index,iatno,iadef(ind),inpcrdtyp,
                call getint('Valence',7,nval(iatnoadd),1,8,newn,0)
                call defaultbondl(nneig(ia1)+1,newn,iatno(ia1),iatnoadd,r12)
                call getreal('R1-X distance',13,r12,rij,1,0)
                call addatom(nndef,c(1,ia4),c(1,ia3),c(1,ia2),c(1,ia1),
    c             Make room
                call blankout(line(index(n)),iresncol1,iresncol2)
                call writeint(line(index(n)),icolw,iresno(n),lenw)
                call rightadjustline(line(index(n)),iresncol1,iresncol2)
                  call blankout(line(index(n)),iseqncol1,iseqncol2)
                  call writeint(line(index(n)),icolw,lastseq+1,lenw)
                  call rightadjustline(line(index(n)),iseqncol1,iseqncol2)
    c         Add amide hydrogens
              call getint('Index of the first atom of the new bond',39,
              call listatom(line,index,iatno,ia1,inpcrdtyp,iofull,
              call getint('Index of the second atom of the new bond',40,
              call listatom(line,index,iatno,ia2,inpcrdtyp,iofull,
    c         Specify properties of the new atom
                call getname(pflab,len,question,lq,ipotcol2-ipotcol1+1,'',0,
                  chargech=chargech-qprev
                call getreal('Charge',6,999999.0,q,0,0)
                call putreal(line(index(ia))(iqcol1:iqcol2),iqcol2-iqcol1+1,
                chargech=chargech+q
                charge(ia)=q
                call getname(ATD4typ,len,'Autodock-4 atom type',20,2,'',0,0,
          character*1 chain,ans
          character*3 resnam3
          character*4 atnam4
          character*8 atnam,resnam,resnamn,hnname
          character* 132 line(maxrec)
    c     print *,'RETROINVERSO nslt=',nslt,' iout=',iout
          call quiz(ans,itertyp,' ',' ',0,'terminus type',13,0,5,6,00)
          call indexit(mask,1,nslt,0)
    c     print *,'IAT1,2=',iat1,iat2
            call leftadjust4(atnam(1:4),atnam4)
    c       print *,'IA=',ia,' ATNAM=',atnam(1:lnam)
              call trnsfr(cc,c(1,ia),3)
              call trnsfr(co,c(1,ia),3)
              call trnsfr(cn,c(1,ia),3)
              call trnsfr(ca,c(1,ia),3)
              call trnsfr(ch,c(1,ia),3)
    c         print *,'IRESN=',iresn
    c         New residue
    c         print *,'NRES=',nres
    c         Now make the transformation
    c           Remove extra hydrogens from N (if any)
                  call leftadjust4(atnam(1:4),atnam4)
    c             print *,'IN=',in,' atno=',iatnum(in),' nhf=',nhf
    c           print *,'H clean ndel=',ndel_extra
                  call addatom(2,cg,c(1,ineig(2,infound)),
                    call leftadjust4(atnam(1:4),atnam4)
    c             Generate the amide H
                      ccnk=(cn(k)+ca(k))/2.0
                      co1(k)=2.0*ccnk-cc(k)
                      call leftadjust4(atnam(1:4),atnam4)
                        call trnsfr(co1,c(1,iaa),3)
                  call addatom(2,cg,co1,ca,cn,ch,rnh,0.0,0.0,0.0,1,pi,1,
    c         print *,'NRES=',nres,' IA1,2=',ia1,ia2
    c           Remove extra oxygen from C (if any)
                  call leftadjust4(atnam(1:4),atnam4)
    c             print *,'IN=',in,' name=',atnam(1:4)
    c           Remove ACE terminus atoms (if any)
                  call leftadjust4(atnam(1:4),atnam4)
                    call leftadjust4(atnam(1:4),atnam4)
    c             If Amber Nter was read, search for ca,ch,ch 
                      call leftadjust4(atnam(1:4),atnam4)
                        call trnsfr(cc,c(1,i),3)
                        call trnsfr(co,c(1,i),3)
                        call trnsfr(ca,c(1,i),3)
    c       print *,'CA=',ca
    c       print *,'CO->H=',co
    c       print *,'CC->N=',cc
    c       print *,'RESNAM3=',resnam3
    c       print *,'RESNAM3=',resnam3,' IAMBER=',iamber
    c           print *,'IA1,2=',ia1,ia2,' IAMBER=',iamber
    c           print *,'O clean ndel=',ndel_extra
    c             Zwitterion - generate H1,H2,H3
    c             print *,'DANG123=',dang1,dang2,dang3
                  call addatom(1,cg,cn,ca,cc,ch23,rnh,ang,dang2,0.0,1,pi,1,
                  call addatom(1,cg,cn,ca,cc,ch23,rnh,ang,dang3,0.0,1,pi,1,
    c             Blocked - generate ACE as separate residue (Amber)
                  call addatom(2,cg,ca,co,cc,cco,rcn,0.0,0.0,57.0,1,pi,1,
                  call addatom(1,cg,ca,cc,cco,oco,rco,120.0,90.0,0.0,1,
                  call addatom(2,cg,oco,cc,cco,ccc,rcc,0.0,0.0,0.0,1,
                  call addatom(1,cg,cn,cco,ccc,ch1,rch,104.0,0.0,0.0,1,
                  call addatom(1,cg,cn,cco,ccc,ch2,rch,104.0,120.0,0.0,1,
                  call addatom(1,cg,cn,cco,ccc,ch3,rch,104.0,-120.0,0.0,1,
    c             Blocked - generate ACE as same residue
    c       print *,'CA=',ca
    c       print *,'CC=',cc
    c       print *,'CO=',co
    c       print *,'RESNAM3=',resnam3
                  call addatom(2,cg,ca,co,cc,cco,rcn,0.0,0.0,57.0,1,pi,1,
                  call addatom(1,cg,ca,cc,cco,oco,rco,120.0,90.0,0.0,1,
                  call addatom(2,cg,oco,cc,cco,ccc,rcc,0.0,0.0,0.0,1,
                  call addatom(1,cg,cn,cco,ccc,ch1,rch,104.0,0.0,0.0,1,
                  call addatom(1,cg,cn,cco,ccc,ch2,rch,104.0,120.0,0.0,1,
                  call addatom(1,cg,cn,cco,ccc,ch3,rch,104.0,-120.0,0.0,1,
    c           Remove N terminal H and/or NME
                  call leftadjust4(atnam(1:4),atnam4)
    c         write (6,8932) ia1,ia2,nres,iresn,iamber
    c8932     format(' IA1,2=',2i5,' NRES,IRESN=',2i5,' IAMBER=',i2)
    c           write (6,9861) nres,ia1,ia2,resnam(1:3),ireso,iresn,
    c    -        icfound,infound,iofound,ihfound,iafound
    c9861       format(' nres=',i4,' ia1,2=',2i6,' resnam=',a,' ireso,n=',
    c    -        2i4,/,' icfound,infound,iofound,ihfound,iafound=',5i6)
    c             Special procedure for proline
    c             print *,'PRO ia1,ia2=',ia1,ia2
                    call leftadjust4(atnam(1:4),atnam4)
                      call trnsfr(ch,c(1,ihfound),3)
                  call trnsfr(cn,c(1,icfound),3)
                  call trnsfr(ca,c(1,iafound),3)
                  call trnsfr(cb,c(1,ibfound),3)
                  call trnsfr(cd,c(1,iofound),3)
    c              write (6,9866) icfound,iafound,ibfound,iofound,infound
    c9866          format(' icfound,iafound,ibfound,iofound,infound=',5i6)
    c             write (6,*) 'cn=',cn
    c             write (6,*) 'ca=',ca
    c             write (6,*) 'cb=',cb
    c             write (6,*) 'cd=',cd
    c             print *,'phirad=',phirad,' ',phirad
    c             Generate new CG
                    call addatom(1,cd,cn,ca,cb,cg,1.0,105.0,phirad,0.0,0,pi,
                      cm_na(k)=(c(k,infound)+ca(k))/2.0
    c               write (6,*) 'cm_na=',cm_na
                      cm_na(k)=(cn(k)+ca(k))/2.0
                      cm_bd(k)=(cb(k)+cd(k))/2.0
    c               write (6,*) 'cm_na=',cm_na
    c               write (6,*) 'cm_bd=',cm_bd
                    call arrdiff(cm_bd,cm_na,e,3)
    c               write (6,*) 'e=',e
    c               write (6,*) 'dcg=',dcg,' enorm=',enorm
                      cg(k)=cm_na(k)+dcg*e(k)/enorm
    c             write (6,*) 'cg=',cg
    c             Generate new hydrogens
                  call addatom(2,e,ca,cg,cb,rh,1.0,0.0,0.0,57.,0,pi,1,ifail)
                  call addatom(2,e,ca,cg,cb,rh,1.0,0.0,0.,-57.,0,pi,1,ifail)
                  call addatom(2,e,cb,cd,cg,rh,1.0,0.0,0.0,57.,0,pi,1,ifail)
                  call addatom(2,e,cb,cd,cg,rh,1.0,0.0,0.,-57.,0,pi,1,ifail)
                  call addatom(2,e,cn,cg,cd,rh,1.0,0.0,0.0,57.,0,pi,1,ifail)
                  call addatom(2,e,cn,cg,cd,rh,1.0,0.0,0.,-57.,0,pi,1,ifail)
    c             Write new backbone
                  call rescale_bl(c(1,icfound),c(1,iofound),chn,rco,rnh)
                  call rescale_bl(c(1,infound),ch,con,rnh,rco)
    c         print *,'IA1,IA2,IA,NRES=',ia1,ia2,ia,nres
    c       Zwitterion - generate O2
            call addatom(2,cg,con,c(1,iaf),c(1,inf),e,rco,0.0,0.0,
    c       Blocked - generate NME as separate residue
    c       print *,'CA=',ca
    c       print *,'CN=',cn
    c       print *,'CH=',ch
            call addatom(2,cg,ca,con,cn,cnn,rcn,0.0,0.0,0.0,1,pi,1,ifail)
            call addatom(1,cg,ca,cn,cnn,cnh,rnh,104.0,-120.0,0.0,1,
            call addatom(2,cg,cn,cnh,cnn,ccc,rcn,0.0,0.0,0.0,1,pi,1,ifail)
            call addatom(1,cg,cn,cnn,ccc,ch1,rch,104.180,0.0,0.0,1,pi,1,
            call addatom(1,cg,cn,cnn,ccc,ch2,rch,104.0,120.0,0.0,1,pi,1,
            call addatom(1,cg,cn,cnn,ccc,ch3,rch,104.0,-120.0,0.0,1,pi,1,
    c       Blocked - generate NME as same residue
            call addatom(2,cg,ca,con,cn,cnn,rcn,0.0,0.0,0.0,1,pi,1,ifail)
            call addatom(1,cg,ca,cn,cnn,cnh,rnh,104.0,-120.0,0.0,1,
            call addatom(2,cg,cn,cnh,cnn,ccc,rcn,0.0,0.0,0.0,1,pi,1,ifail)
            call addatom(1,cg,cn,cnn,ccc,ch1,rch,104.0,0.0,0.0,1,pi,1,ifail)
            call addatom(1,cg,cn,cnn,ccc,ch2,rch,104.0,120.0,0.0,1,pi,1,
            call addatom(1,cg,cn,cnn,ccc,ch3,rch,104.0,-120.0,0.0,1,pi,1,
            cn(k)=c1(k)+e(k)*rn
    c     Establish default bondlength based on atomic numbers and valences
    c     To be extended
    c       Carbon
    c       if (nngi1 .eq. 4) then
    c       else if (nngi1 .eq. 3) then
    c         if (iatn2 .eq. 1) r12=1.03
    c         if (iatn2 .eq. 6) r12=1.03
    c         if (iatn2 .eq. 1) r12=1.03
    c         if (iatn2 .eq. 6) r12=1.03
    c       end if
    c       Nitrogen
    c       if (nngi1 .ge. 3) then
    c       end if
    c       Oxigen
    c         if (iatn2 .eq. 7) r12=1.14
    c         if (iatn2 .eq. 8) r12=1.03
    c         if (iatn2 .eq. 7) r12=1.03
    c         if (iatn2 .eq. 8) r12=1.03
    c       Phosphorus
    c       if (iatn2 .eq. 1) r12=1.03
    c       if (iatn2 .eq. 6) r12=1.03
    c       if (iatn2 .eq. 7) r12=1.03
    c       if (iatn2 .eq. 8) r12=1.03
    c       Sulphur
    c       if (iatn2 .eq. 7) r12=1.03
    c       if (iatn2 .eq. 8) r12=1.03
          character*2 lev
          character* 132 line(maxrec)
          character*16 lab1,lab2
    c     print *,'FINDROOT ia=',ia,' inpcrdtyp=',inpcrdtyp,' ix(ia)=',index(ia)
    c       print *,'ia,ja,ia1,iatno(ia1)=',ia,ja,ia1,iatno(ia1)
              call listatom(line,index,iatno,ia1,inpcrdtyp,iofull,lab1,16,n,
          character* 132 line(maxrec),ll
          character*(*) lab
          character*4 pflab
          character*12 qlab
          call setcol(inpcrdtyp,ncol,idcol,ialtcol,iinscol,
    c     Add one atom with known distance, angle and torsion
          character*13 lab
    c     print *,'ADDATOM ityp=',ityp
    c       Seqential definition
            call arrdiff(r1,r2,c12,3)
            call arrdiff(r3,r2,c23,3)
    c       print *,'ADDATOM r1=',r1
    c       print *,'ADDATOM r2=',r2
    c       print *,'ADDATOM r3=',r3
                c12(k)=c12(k)/rnorm12
                c23(k)=c23(k)/rnorm23
              call vprd(c12,c23,e3)
              call vprd(e3,c12,e2)
    c         Atoms r3, r2, r1, and x are colinear
    c       print *,'ADDATOM x=',x
            call trnsfr(p(1,1),r3,3)
            call trnsfr(p(1,2),r2,3)
            call trnsfr(p(1,3),r1,3)
            call trnsfr(p(1,4),x,3)
    c       print *,'2,3,4 angle=',angleijk(p,4,2,3,4,6)*180.0/3.1415
    c       print *,'1,2,3,4 dih angle=',
    c    -    dihangl(p,1,2,3,4,0,100000)*180.0/3.1415
    c       Bisector definition
    c       print *,'BEND=',bend,' AIJK,TIJKL=',aijk,tijkl
    c       print *,'R1=',r1
    c       print *,'R2=',r2
    c       print *,'R3=',r3
    c       print *,'RIJ=',rij
            call arrdiff(r1,r2,c12,3)
            call arrdiff(r1,r3,c13,3)
              c12(k)=c12(k)/rnorm12
              c13(k)=c13(k)/rnorm13
            colinear=scprod(c12,c13)
            call arrsum(c12,c13,bs,3)
              call arrdiff(r2,r3,c23,3)
              call vprd(c23,bs,e3)
    c       print *,'X =',x
    c       Trisector definition
            call arrdiff(r1,r2,c12,3)
            call arrdiff(r1,r3,c13,3)
            call arrdiff(r1,r4,c14,3)
              c12(k)=c12(k)/rnorm12
              c13(k)=c13(k)/rnorm13
              c14(k)=c14(k)/rnorm14
    c     Extract parts of the system
          character*1 asterisk,ans,lastans
          character*4 chainid,chainidn,atnam,atnamadj,resnam,resnamadj,
          character* 132 line(maxrec)
    c      write (6,7877) isegcol1,isegcol2,iseqncol1,iseqncol2,
    c     -  inamcol1,inamcol2,irescol1,irescol2
    c7877  format(' isegcol12=',2i4,' iseqncol12=',2i4,' inamcol12=',2i4,
    c     -  ' irescol12=',2i4)
            call readint(line(index(nslt)),iseqncol1,iseqncol2,iseqno,1,1,
          call zeroiti(indexdel,0,n)
            call quiz(ans,ians,' ',' ',0,'selecting option',16,0,5,6,0)
                call blankout(chainid,1,nsegcol)
                  call getname(chainid,len,'Chain ID to keep',16,4,'',0,0,0,
                  call getname(chainid,len,'Chain ID to drop',16,4,'',0,0,0,
                call askstop(1)
                call getrange(ifst,999999,ilst,999999,i,0,'atom to keep',
              call getrange(ifst,999999,ilst,999999,i,0,'atom to drop',12,
                call askstop(1)
                call getrange(ifst,999999,ilst,999999,i,0,
              call getrange(ifst,999999,ilst,999999,i,0,
    c         Read atom name list to keep
              call getnamelist(namesel,4,nnamesel,'Atoms names to use',18,
    c         Read residue name list to keep
              call getnamelist(namesel,nrescol,nnamesel,
    c         Drop all solvents
                call askstop(-1)
                call askstop(-1)
                chainidn=line(index(ia))(isegcol1:isegcol2)
                chainidn='    '
                call readint(line(index(ia)),iseqncol1,iseqncol2,iseqno,1,1,
                call askyn('Do you want to ignore the sequence number read',
              call leftadjust4(atnam,atnamadj)
    c          write (77,7778) ia,iseqno,iresno,atnam,atnamadj,
    c     -      chainid(1:nsegcol),chainidn(1:nsegcol)
    c7778      format(' ia,iseqno,iresno=',3i4,' atnam=',a,' atnamadj=',a,
    c     -      ' chainid,n=',a,1x,a)
    c           if (atnamadj(1:3) .ne. 'CA ' .and. atnamadj(1:2) .ne. 'C '
    c    -          .and. atnamadj(1:2) .ne. 'O '
    c    -          .and. atnamadj(1:2) .ne. 'N '
    c    -          .and. atnamadj(1:2) .ne. 'H '
    c    -          .and. atnamadj(1:3) .ne. 'HN ') idrop=1
    c            write (77,9781) ia,atnam,atnamadj(1:3),idrop
    c9781        format(i6,' atnam=',a,' atnamadj=',a,' idrop=',i2)
                call leftadjustn(resnam,resnamadj,nrescol)
                    call askyn('Do you want to keep it',22,1,-1,ikeep,0,0)
    c         if (idrop .gt. 0 .and.
    c    -       line(index(ia))(idcol:idcol) .ne. asterisk) then
    c           nrecdel=nrecdel+1
    c           line(index(ia))(idcol:idcol)=asterisk
    c         end if
    c         Write back charges into line
                call blankout(line(index(ia)),iqcol1,iqcol2)
          character*(*) question
          character*80 line
    c     Input a range
    c     print *,'GETRANGE ifst,ifstdef,ilst,ilstdef=',
    c    - ifst,ifstdef,ilst,ilstdef
    c     print *,'GETRANGE maxval=',maxval
          call getint(line,6+lquestion,ifstdef,1,maxval,ifst,ihelp)
          call getint(line,6+lquestion,ilstdef,1,maxval,ilst,ihelp)
          character* 132 line(maxrec)
    c     Change atomnumbers in bond list
    c         if (iold .ne. 0) write (6,888) ia,index(ia),ib,ic1,ic2,iold
    c888      format(' ia,index(ia),ib=',3i3,' ic1,2=',2i3,' iold=',i3)
    c           print *,'isort,indexo(iold),indexn(iold)=',
    c    -         isort,indexo(iold),indexn(iold)
    c     Find intg in the sorted list
    c#    MMC routine 127 lstmod: 01/08/88
    c*****Merge two sets for the sorting
          character* 132 line(maxrec)
    c     i,j: counters for the first and second part, resp.
    c     k: counter for the merged array
    c     print *,'m1,m2,n1,n2,n=',m1,m2,n1,n2,n
    c       Pick value from first or second part, depending on the comparison
    c         Compare residue numbers since segment id's were identical
              call readint(line(index(i)),iresncol1,iresncol2,numi,2,1,
              call readint(line(index(j)),iresncol1,iresncol2,numj,2,1,
    c         Leftover in the first part
    c         Lefotver in the second part
    c     Move back merged data into index
    c*****Sort the indexx and value array in the order of increasing value(i)
    c     print *,'MRGSRT n,maxt=',n,maxt
          call mrglimtst(iout,n,maxt,ireturn)
          call indexit(ifa,1,n,0)
          call indexit(ila,1,n,0)
    c     Merge pairs of intervals
            call mergelst(indexx,value,ifa(l),ila(l),ifa(l+1),ila(l+1),
    c       Take care of last (odd) interval
            call mergelst(indexx,value,ifa(nnpair),ila(nnpair),ifa(nn),
    c#    MMC routine 263 lstmod: 06/23/97
    c*****Merge two sets for the sorting
    c     i,j: counters for the first and second part, resp.
    c     k: counter for the merged array
    c       Pick value from first or second part, depending on value
    c       Leftover in the first part
    c       Leftover in the second part
    c     Move back merged data into index, value
    c#    MMC routine 262 lstmod: 10/17/97
    c*****Sort the indexx and ivalue array in the order of increasing ivalue(i)
          call mrglimtst(iout,n,maxt,ireturn)
    c     Individual intervals consist of single elements
    c     print *,'MRGST n,maxt=',n,maxt
          call indexit(ifa,1,n,0)
          call indexit(ila,1,n,0)
    c     Merge pairs of intervals
            call mergelsti(indexx,ivalue,ifa(l),ila(l),ifa(l+1),ila(l+1),
    c       Take care of last (odd) interval
            call mergelsti(indexx,ivalue,ifa(nnpair),ila(nnpair),ifa(nn),
    c#    MMC routine 263 lstmod: 06/23/97
    c*****Merge two sets for the sorting
    c     i,j: counters for the first and second part, resp.
    c     k: counter for the merged array
    c       Pick ivalue from first or second part, depending on ivalue
    c       Leftover in the first part
    c       Leftover in the second part
    c     Move back merged data into index, value
    c     Individual intervals consist of single elements
    c     Find out if any of the stored datapaths has a pdb_nam.dat
          character*100 datapath,datapaths
          common /environment/ npaths,ldatapath,ldatapaths(5),
    c     Try stored paths
            call trypath(datapaths(ip),ldatapaths(ip),datapath,ldatapath)
            call blankout(datapaths(5),1,100)
            call getenv('PWD',datapaths(5))
            call lastchar(datapaths(5),ldatapathtry,100)
            call trypath(datapaths(5),ldatapathtry,datapath,ldatapath)
            call trypath(datapaths(5),ldatapathtry,datapath,ldatapath)
              call askyn('Do you want to try a different path',35,1,1,itry,
          character*100 datapathtry,datapath
          character*200 filename
          call openfile(88,0,' ',1,'old',filename,namlen,notfound,2,1,0,1,
            close (88)
          character* 132 line(maxrec)
          character*8 resnam,resnames(maxrsd),resnamslv
          character*4 segnam
          character*5 crdext
          common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
          common /solvres/ iaskdiffres,idiffres
          common /logging/ logfile,ipredict
    c     Make sure residue numbers and atom numbers are consecutive
          call setcol(inptyp,ncol,idcol,ialtcol,iinscol,
    c     print *,'RESEQ nslt,naslv,nconfig=',nslt,naslv,nconfig
    c     print *,'reseq  n=',n,' nrescol,nresncol,nsegcol=',
    c    -  nrescol,nresncol,nsegcol
    c     print *,irescol1,irescol2,iresncol1,iresncol2,
    c    -  isegcol1,isegcol2
    c     print *,' iseqncol1,2=',iseqncol1,iseqncol2
    c     print *,'iresidcol1,iresidcol2=',iresidcol1,iresidcol2
            call askyn('Do you want to adjust atom and residue numbers',
              call getint('Atom number of the first atom',29,
              call getint('Residue number of the first residue',35,
    c         print *,'First atom and residue numbers=',iatfirst,iresfirst
                call getint('Residue id of the first residue',31,
    c           Allow for consecutive residue numbers over segments
                call askyn(
                call askyn(
    c     print *,'resnam,resnum,segnam,n=',resnam,resnum,segnam,n
    c     print *,'C nsegcol=',nsegcol
    c           Second solvent is the same residue as the first
                call askyn(
    c           New segment was found
                  call getint('Residue number of the first residue',35,
                  call getint('Residue id of the first residue',31,
    c       Sequence number (if used)
              call readint(line(index(ia)),iseqncol1,iseqncol2,iaorg,1,1,
    c       Check for change in residue number/id
            call readint(line(index(ia)),iresncol1,iresncol2,iiresorg,2,1,
    c         Residue id
              call readint(line(index(ia)),iresidcol1,iresidcol2,iiresorg,2,
    c     Add '3' to 'TIP'
    c     print *,'numres=',numres
          character* 132 line(maxrec)
          character*4 atnamo,atnamn
          character*8 resnam
    c     'Regularize' PDB atomnames
    c     'De-regularize' PDB atomnames
            call regularpdb(atnamo,atnamn,itofrom)
          character*132 line(maxrec)
          character*1 altcol(maxrec),inscol(maxrec)
          character*8 atnames(maxrec)
          character*1 asterisk
          character*5 crdext
          common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
    c     print *,'ATOMDEL nrecdel,n,idcol,nconfig=',nrecdel,n,idcol,nconfig
    c       write (77,*) 'ATOMDEL i=',i,' AC:',line(index(i))(idcol:idcol)
    c       Update index to eliminate the dropped molecules
              charge(i-ndel)=charge(i)
              call trnsfr(c(1,i-ndel),c(1,i),3)
    c     print *,'atomdel n,inpcrdtyp=',n,inpcrdtyp
    c     Establish atom ranges for the residues reflecting the deletions
    c     Read in atom and residue name conversion tables, set up limits, etc
          character*200 filename,line
          character*8 namnam
          character*9 recnam
          common /convspec/ incon,ioutcon,ideoxy,namnam(10)
          character*3 ires0
          character*4 anamcontab,rnamcontab,atnam
          common /convdat/ nanamcon,nrnamcon,nres0,ifst(100),ilst(100),
          call openfile(iuconv,0,'namcon',6,'old',filename,namlen,
    c     Read in conversion information
    c     nconcol: number of colums containing conversion data
    c     ncon: number of name conversions (PDB can have several columns)
    c     namename: Names of the conventions
          call chosename(namnam,ncon,'input',5,incon)
          call chosename(namnam,ncon,'output',6,ioutcon)
    c       read (iuconv,*,err=9921,end=992)
    c       read (iuconv,*,end=992)
    c    -    (rnamcontab(j,i),j=1,nconcol)
    c        write (77,2200) i,irec,(anamcontab(irec+1,j),j=1,nconcol)
    c2200    format(' i,irec=',2i4,' anam=',12(1x,a4))
    c        write (77,7711) i,irec,nconcol,ndiff,anamcontab(irec+1,1),
    c     -   anamcontab(irec+1,2)
    c7711    format(' i,irec,nconcol,ndiff=',4i5,2a6)
    c         New residue found
          close (iuconv)
    c      do i=1,nres0
    c        write (99,7123) i,ires0(i),ifst(i),ilst(i)
    c7123    format(i3,1x,a3,3i5)
    c      end do
    c      do i=1,nanamcon
    c        write (77,7755) i,(anamcontab(i,j),j=1,nconcol)
    c7755    format(i5,10(2x,'*',a4,'*'))
    c      end do
          character*(*) namelist(nnames),lab
          call getint('Convention number',17,999999,1,nnames,iname,00)
    c     Perform an atom and residue name conversion
          character* 132 line
          character*1 ans1
          character*4 atnamin,atnamout,an0,an1,resnamin,resnamout,rn0
          character*8 namnam
          common /convspec/ incon,ioutcon,ideoxy,namnam(10)
          character*3 ires0
          character*4 anamcontab,rnamcontab
          common /convdat/ nanamcon,nrnamcon,nres0,ifst(100),ilst(100),
    c     write (77,*) 'Aresnamin=',resnamin,' atnamin=',atnamin
          call leftadjust4(resnamin,resnamin)
          call leftadjust4(atnamin,an1)
            call regularpdb(an1,atnamin,1)
    c     write (77,*) ' an1=',an1,' atnamin=',atnamin
    c     print *, 'an1,atnamin=',an1,'|',atnamin
    c     Find first the residue name in the resname table
    c           write (77,*) 'Match i,j=',i,j,' rno=',resnamout
    c     Find generic residue name
    c         write (77,*) 'Generic resname=',ires0(i)
    c     Find the atom conversion - first check residue-specific conversions
    c           write (77,*)
    c    -        'i=',i,' atnamin=',atnamin,' atab=',anamcontab(i,j)
    c             write (77,7712) atnamin,resnamin,atnamout
    c7712          format(' FOUND: an=',a4,' rn=',a4,' ao=',a4)
    c     Not found among the specifics  - check the residue independent list
    c           write (77,*) 'Generic found atnamout=',atnamout,' i,j=',i,j
    c     Take care of oxy/deoxy nucleic acids
    c             Ask if oxy or deoxy NA
                  call getname(ans1,len,
    c       Gromacs output - check for residue name change
              call findname('HD2 ',anamcontab(1,ioutcon+1),ifst(iresconv0),
              call findname('HG  ',anamcontab(1,ioutcon+1),ifst(iresconv0),
              call findname('HE2 ',anamcontab(1,ioutcon+1),ifst(iresconv0),
              call findname('HZ3 ',anamcontab(1,ioutcon+1),ifst(iresconv0),
    c     Check for changes
    c         Undo the regularization
              call leftadjust4(an0,an1)
    c      write (77,7723) rn0,an0,resnamout,atnamout
    c7723  format(' rn0,an0=',a,1x,a,' ,resnamout,atnamout=',a,1x,a)
          character*4 inp,reg,regg,inpl
          call leftadjust4(inp,inpl)
    c     write (06,*) 'REGPDB inp=',inp,' inpl=',inpl,' itofrom=',itofrom
    c       To PDB
    c           Two-character H name
    c           Three-character H name
    c           Four-character H name
    c             HXXd -> dHXX
    c             HddX -> dHXd
    c             HdXX -> dHXX
    c             HXdX -> dHXX
    c       From PDB
    c       Just move the number away from the first position, if any
    c           Two-character H name
    c           Three-character H name
    c           Four-character H name
    c           Permute cyclicly
    c     write (06,*) 'REGPDB reg=',reg
          character*132 line(maxrec)
          character*1 lcprev1,lcprev2
          character*5 atomnam,prev1,prev2
    c     Look for atomnames X2,X3 without X1
            call lastchar(atomnam,lc,nnamcol)
    c        write (77,2299) ia,atomnam,lc,lcprev1,lcprev2,prev1,prev2
    c2299    format(i5,' anam=',a,' lc=',i1,' lcprev1,2=',a,1x,a,
    c     -    ' prev1,2=',a,1x,a)
    c           X2,X3 without X1 found
                  call askyn('Do you want to shift H*2-H*3 to H*1-H*2',39,
    c           Make the shift
          character*(*) prompt,filename
          character*3 mode
          character*100 datapath,datapaths
          common /environment/ npaths,ldatapath,ldatapaths(5),
          common /logging/ logfile,ipredict
          character*11 formnam
          character*80 promptq
    c     iswitch=1: Switch to terminal input if namlen=0
    c     ib: 1/2: formatted/unformatted
    c     nosys=1: don't try to look into datapath
    c     ioverall=1: no matter what, overwrite without asking existing file
    c     idatapath=1: when an 'old' file is not found, try it in the datapath dir.
    c     idatapath=2: when an 'old' file is not found, just returnt nofund=1
    c     idatapath=3: when an 'old' file is not found, ask for a new file name
    c     print *,'OPENFILE ipredict=',ipredict,' fn=',filename(1:namlen)
    c     print *,'OPENFILE ib,f,=',ib,formnam(ib),' fn=',filename(1:namlen)
    c     print *,'OPENFILE namlen=',namlen
            call getname(filename,namlen,promptq,lprompt,200,'',0,0,0,0)
    c       print *,'Name read:',filename(1:namlen)
            call askyn('If the file exists, do you want to overwrite it',47,
    c         Just return notfound=1
    c           Try datapath
                call askyn(
                  close (iunit,status='delete')
          character*(*) dirname
          character*200 filename
          character*(*) oldname,newname,newext
          character*80 question,oldext
          character*200 tempname
    c       Input file had no '.' in it - just add extension
    c       Replace the extension
              call askyn(question,lq,1,1,ians,00,0)
    c       Insert extension between file name root and old extension
    c     idcdtyp: 0: not a DCD file; 1: coordinate DCD file; 2: velocity DCD file
          character*200 filename
          character*4 header
          call openfile(97,0,' ',1,'old',filename,namlen,notfound,0,2,1,1,0)
    c     itrajtyp: 0: not aan Amber trajectory file
    c               1: likely to be an Amber tarjectory file
          character*200 filename
          character* 132 line
          call openfile(97,0,' ',1,'old',filename,namlen,notfound,0,1,1,1,0)
          call blankout(line,1,132)
          call lastchar(line,ifc,132)
    c    -  filename,namlen,notfound,idatapath,ib,nosys,noecho,ioverall)
          character*8 resnames(maxrsd)
          character*1 resnames1(maxrsd)
          character*200 inpfile,outfile
          character*4 segnames(maxrsd) 
          character*1 ans,ch
          character*4 s4
          character*132 l132(1)
          character*200 line
          character*1 abc,digits,hexdigits
          common /charactersets/ ihex(25),abc(62),digits(14),hexdigits(25)
          call openfile(30,0,'input sequence',14,'old',inpfile,linpfile,
          call quiz(ans,inpstyp,' ',' ',0,
          call openfile(31,0,'output sequence',15,'new',outfile,loutfile,
          call quiz(ans,iostyp,' ',' ',0,
    c       Charmm sequence input
              call blankout(line,1,80)
              call lastchar(line,lc,80)
    c           Header read
                  call blankout(line,1,80)
                  call lastchar(line,lc,80)
                    call changeprot(resnames(nres)(1:3),resnames1(nres),2)
                call blankout(line,1,80)
    c             print *,'GENERATE * SETUP record is missing'
                  call nextblank(line,ic,80)
    c       PDB SSEQRES input
            ch=' '
              call blankout(line,1,80)
                call lastchar(line,lc,80)
                ch=line(12:12)
    c           print *,'L12=',line(12:12),' CH=',ch,' NS=',nseg,' NR=',nres
                  call changeprot(resnames(nres)(1:3),resnames1(nres),2)
    c       >Title + <1>-char list (FASTA) or PIR
            call blankout(line,1,80)
              call blankout(line,1,80)
                call lastchar(line,lc,80)
                  call changeprot(resnames(nres)(1:3),resnames1(nres),1)
    c       CIF
              call blankout(line,1,200)
              call lastchar(line,lc,200)
              call findnextchar('p',line,ic,200)
                  ch=abc(nseg+1)
    c               Sequence in the same line
                      call nextblank(line,ic,200)
                    call nextchar(line,ic,200)
                      call changeprot(resnames(nres)(1:3),resnames1(nres),1)
                      call blankout(line,1,200)
                        call lastchar(line,lc,200)
                          call changeprot(resnames(nres-incr+ic)(1:3),
    c       GCG
              call blankout(line,1,200)
              call lastchar(line,lc,200)
                call nextchar(line,ic,200)
                call nextblank(line,ic,200)
                call nextchar(line,ic,200)
                call nextblank(line,ic,200)
                  call changeprot(resnames(i)(1:3),resnames1(i),1)
    c           print *,'NSEG=',nseg,' NRES=',nres,' S4=',s4
    c       print *,'ISEG=',iseg,' INPSTYP=',inpstyp,' IOSTYP=',iostyp
    c         Charmm sequence format
    c         PDB SSEQRES iout
    c         >Title + <1>-char list (FASTA)
    c         PIR
    c         GCG
              call gcgwrite(0,l132,0,outfile,loutfile,31,
          close(30)
          close(31)
    c*****Generates a residue list in a variety of formats
          character*200 outfile,seqfile
          character* 132 line(maxrec)
          character*80 title
          character*1 ans,resnames1(maxrsd)
          character*4 segid
          character*8 resnames(maxrsd)
          character*1 aanames1
          character*2 mmodtoamb
          character*3 aanames3
          character*80 seqid
          common /atnamcon/ mmodtoamb(100),aanames1(58),aanames3(58),
          common /columnlim/ incol(19),iidcol(19),iialtcol(19),iiinscol(19),
          character*1 abc,digits,hexdigits
          common /charactersets/ ihex(25),abc(62),digits(14),hexdigits(25)
          character*5 crdext
          common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
          character*6 typnam(5)
          character*30 question
    c     print *,'WRITESEQ numres,numslv=',numres,numslv
    c       Run incidental to sequence generation - use default format w/o/ asking
            call askyn('Do you want a sequence list',27,1,-1,iseqls,0,0)
    c       Ask for format
            call quiz(ans,iostyp,' ',' ',0,
          call changeext(outfile,seqfile,namleno,namlens,'seq',3,0,0)
          call openfile(30,0,' ',1,'new',seqfile,namlens,notfnd,0,1,1,0,0)
    c     Write sequence input (by segments)
              call nextchar(title,ifc,80)
    c         Convert to 1-character ID-s
                call changeprot(resnames(ir),resnames1(ir),2)
    c           PIR
                call getname(seqid,namlen,question,25,80,abc(nseg),1,0,0,0)
                call lastchar(title,lentit,80)
                call gcgwrite(inpcrdtyp,line,index(1),seqfile,namlens,30,
          close (30)
          character*1 seq(lenseq),type
          character*4 pdbid
          character* 132 line(maxrec)
          character*200 inpfile
          character*40 keywords
          character*62 textinp
          character*5 crdext
          common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
    c     print *,'GCGWRITE lenseq=',lenseq,' ifile=',ifile
    c     Determine if protein or nucleic acid
    c         Amino acid residue found
    c         Translate header into keywords and pdb id
            call getname(keywords,len,'Keywords (max 40 chars)',23,40,'',0,
            call getname(pdbid,len,
    c     Calculate the checksum
    c         Get compound name as definition
            call getname(textinp,len,'Compound name (for DEFINITION key)',
          common /pbcrotmat/ torot_ac(3,3),torot_ca(3,3),tofac_ac,tofac_ca
          common /rotmat/ matrot0(4,4),matrot(4,4),nomat0
          character*1 xyz
          common /axislab/ xyz(3)
          character*1 ans,ansedge,ansvertex
    c       None
    c       Cubic
    c       Rectangular
    c       Face-centered cubic
    c       Face-centered cubic
    c       Hexagonal prizms
            call quiz(ansvertex,ixyzhex(2),' ',
    c         Skewed hexagonal prizm
    c       Truncated Octahedron,
    c         Charmm convention (x normal to square)
    c         Truncated Octahedron, Amber/NAMD convention (X normal to hexagon)
    c       Hexagonal close packed
    c       Octahedral
    c       Image-cell file input
    c       Spehere
          character*1 xyz
          common /axislab/ xyz(3)
          common /numbers/ sq3,sq3inv,sq3p2,sq2p3
    c       Will be read with the images
            call getreal('Edge length (A)',15,999999.0,edge(1),1,0)
              call getreal('Edge length in the '//xyz(k)//' direction (A)',
            call getreal('Edge parameter of the FCC cell (A)',34,999999.0,
    c       Hexagonal prism (regular)
            call getreal('Length of the prism (A)',23,999999.0,edge(1),1,0)
            call getreal('Edge of the hexagon (A)',23,999999.0,edge(2),1,0)
    c       Skewed hexagonal prism
            call getreal('Length of the prism (A)',23,999999.0,edge(1),1,0)
            call getreal('Cell length (a) along Cartesian axis (A)',40,
            call getreal(
    c       Truncated octahedron, Charmm convention (x axis to a square face)
            call getreal('Charmm periodic cell X coordinate (A)',37,
    c       Truncated octahedron, Amber/NAMD convention (x axis to a hexagon face)
            call getreal('Amber/NAMD periodic cell X coordinate (A)',41,
    c       Hexagonal close packing
            call getreal('Inscribed sphere diameter (A)',29,999999.0,
            call getreal('Edge of the rhomboid (A)',24,999999.0,
            call getreal('Radius of the sphere',20,999999.0,edge(1),1,0)
          character*1 xyz
          common /axislab/ xyz(3)
          character*1 abc,digits,hexdigits
          common /charactersets/ ihex(25),abc(62),digits(14),hexdigits(25)
          character*1 ans
          character*80 line
          character*200 imgfile
    c     print *,'READIMG cell=',scale
          call quiz(ans,ians,' ','image',5,'PBC cell type',13,0,5,6,0)
          call openfile(30,0,'image',5,'old',imgfile,lenfilename,notfnd,
    c       Charmm
              call getreal(
            call zeroit(cell,3)
    c         write (6,888) line
    c           write (6,1003) (ic,(cell(k,ic),k=1,3),ic=1,ncell)
                    cell(k,ncell)=-scale(k)*cell(k,ncell)
    c       Simulaid format
              call getreal('Multiplying factor for the '//xyz(k)//' axis',
                cell(k,ncell)=cell(k,ncell)*scale(k)
    c1003  format(' The PBC cell centers:',/,30(i5,3f10.4,/))
          character*(*) quest
            call pbctype(ioppbc,npbc,ixyzhex,nonone)
            call pbcsize(ioppbc,edge,npbc)
              call readimg(cell,ncell)
              call crorgn(edge,edge_gen,ioppbc,3,ncell,cell,cellalt,
            call prtcell(ioppbc,edge,edge_gen,0.0,vol,nw,1)
    c     Calculate the center and radius of the smallest sphere enclosing
    c     all atoms in c
          character*2 iatnm2
          common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
          call zeroit(czero,3)
    c     print *,'COMPACT n,nnh=',n,nnh
          call extension(co,ih,nnh,1,n,cmin,cmax,c0,0,0,v)
    c       print *,'c0,cmax,cmin=',c0(k),cmax(k),cmin(k)
            call cofms(co,com,n,r)
    c     Calculate distances from c0, collect those within 60% of r0
    c     print *,'R02=',r02
    c       print *,'i=',i,' ri2,r02=',ri2,r02
            call askyn('Do you want to keep all atoms instead',37,
    c     print *,'nlist,rorgext=',nlist,rorgext
    c     write (6,2344) (list(i),i=1,6),(r(i),i=1,6)
          call swapl(list,r,nlist,lmax,-1,maxrec)
          call lpshiftc(co,list,nlist,del,c0,np,maxrec)
    c     write (6,2344) (list(i),i=1,6),(r(i),i=1,6)
          call findshift(co,r,list,c0,del,nlist+1,nlist,rlambda,lmax,nzero,
    c     print *,'nzero(1)=',nzero
    c     print *,'lambda=',rlambda
          call transc(c0,del,rlambda)
          call newdist(r,co,list,c0,nlstorg,maxrec)
    c     write (6,2344) (list(i),i=1,6),(r(i),i=1,6)
          call swapl(list,r,nlist,lmax,-1,maxrec)
          call lpshiftc(co,list,nlist,del,c0,np,maxrec)
    c2344 format(' Nlist=',6i2,' r=',6f8.3)
          call findshift(co,r,list,c0,del,nlist+1,nlist-nlistdel,
    c     print *,'nzero(2)=',nzero
    c     print *,'lambda=',rlambda
            c1(k)=co(k,list(nlist+1))-c0(k)
    c       Make these two atoms the diameter
    c     Check the triangle lmax, nlist+1, nlist+2 if it has an obtuse angle
              call swapl(list,r,nlist+ic,lmax,0,maxrec)
              call swapl(list,r,nlist-nlistdel,lmax,0,maxrec)
    c         nlistdel>0 will make the lambda search bypass the last nlistdel ats
              call transc(c0,del,rlambda)
              call newdist(r,co,list,c0,nlstorg,maxrec)
    c         Go back to generate new triangle
          call transc(c0,del,rlambda)
          call newdist(r,co,list,c0,nlstorg,maxrec)
    c     write (6,6666) r(lmax),(r(nlist+ii),ii=1,np)
            call swapl(list,r,nlist,lmax,-1,maxrec)
            call swapl(list,r,nlist,nlist-nlistdel,0,maxrec)
            call swapl(list,r,nlist,lmax,-1,maxrec)
    c     write (6,6666) r(lmax),(r(nlist+ii),ii=1,np)
    c6666  format(' Distances of vertices from c0=',4f8.2)
    c     Find shift preserving the distance from three points
            c1(k)=co(k,list(nlist+2))-co(k,list(nlist+1))
            c2(k)=co(k,list(nlist+3))-co(k,list(nlist+1))
          call vprd(c1,c2,del)
            c1(k)=co(k,list(nlist+1))-c0(k)
          call norm(del,facsgn)
          call findshift(co,r,list,c0,del,nlist+1,nlist-nlistdel,
    c     Check for maximum
    c       First three atoms are in a diagonal plane
    c     Check for center of the sphere being outside the tetrahedron
          call transc(c0,del,rlambda)
          call newdist(r,co,list,c0,nlstorg,maxrec)
    c       iopp(i) will be set to ip if atom (nlist+ip) and c0 are separated
    c       by a tetrahedron face
            call getindex(nlist,ip,iv,i1,i2)
              c1(k)=co(k,list(i1))-co(k,ix)
              c2(k)=co(k,list(i2))-co(k,ix)
              c3(k)=co(k,list(iv))-co(k,ix)
              c4(k)=c0(k)-co(k,ix)
            call vprd(c1,c2,ca)
    c       First drop (temporarily) vertex nlist+iopp(1)
    c       swap vertex ip with nlist+1
            call swapl(list,r,nlist+1,ip,0,maxrec)
    c       rr=sqrt(r(nlist+1))
    c       print *,'R=',rr
    c       write (6,6666) r(lmax),(r(nlist+ii),ii=1,np)
    c       nlistdel>0 will make the lambda search bypass the last nlistdel ats
    c       Go back to check the triangle
    c       Two opposing vertices found
              call getindex(nlist,ip,iv,i1,i2)
              call angles(dist2(co(1,i1),co(1,i2)),r(i1),r(i2),
    c           Acute triangle face found, just drop ip
    c           Swap vertex ip with lmax
                call swapl(list,r,nlist+ip,lmax,0,maxrec)
    c           Move ip to the top of the list for temporary disregard
                call swapl(list,r,nlist,lmax,0,maxrec)
    c           rr=sqrt(r(nlist+1))
    c           print *,'R=',rr
    c           write (6,6666) (r(nlist+ii),ii=1,np)
    c           Go back to generate new tetrahedron
    c       If reached here, neither faces are acute triangles - temporarily
    c       drop both
            call swapl(list,r,nlist+iopp(1),nlist+1,0,maxrec)
            call swapl(list,r,nlist+iopp(2),nlist+2,0,maxrec)
    c     If reached here, center is inside the tetrahedron - done.
          call swapl(list,r,nlist,lmax,-1,maxrec)
    c     Calculation done, perform the shift
          call newdist(r,co,list,c0,nlstorg,maxrec)
    c     write (6,6666) (r(nlist+ii),ii=1,np)
          call shiftmol(co,n,c0,cn,-1.0)
    c     Calculate the minimum radius (as the maximum of the surface atom radii)
          call chksphr(cn,ih,n,nnh,rmin,czero)
    c     Shift the atom lmax to the end of the list,
    c     change nlist by nlincr (-1 or 0)
    c     print *,'nlist,lmax=',nlist,lmax
    c     print *,'SWAPL nlist,lmax=',nlist,lmax
          call swapi4(list(nlist),list(lmax))
    c     Generate linear or planar shift
          call norm(del,1.0)
    c     Calculate the allowed shift of center, rlambda, along del
          coslrl=0.0
            coslrl=coslrl+del(k)*(c(k,list(is))-c0(k))
    c     cosl=cosl/rl
            cosiri=0.0
              cosiri=cosiri+del(k)*(c(k,list(i))-c0(k))
    c       cosi=cosi/ri
    c       if (rl*cosl .ne. ri*cosi) then
    c         rlambdai=(rl2-ri2)/(2.0*(rl*cosl-ri*cosi))
    c        write (6,7777) i,rl2,ri2,rl,ri,cosl,cosi,rlambdai
    c7777    format(i3,' rli2=',2f8.3,' rli=',2f6.2,' cli=',2f8.5,
    c     -    ' l=',f8.3)
              call swapl(list,r,nlist-nzero,lmax,-1,maxrec)
    c     Generate the vertex indices corresponding the (nlist+ip)th list element
            c0(k)=c0(k)+rlambda*del(k)
    c     Check if nothing is outside the sphere
          ca1=(e2+e3-e1)/(2.0*ee2*ee3)
          ca2=(e1+e3-e2)/(2.0*ee1*ee3)
          ca3=(e1+e2-e3)/(2.0*ee1*ee2)
    c     rn is the normal to the plane of a,b, and c
          call vprd(d1,d2,rn)
          call norm(rn,1.0)
    c     r is the radius of the circle going through a,b, and c
          call angdistw(b,a,c,rHB,rb,rac,angabc)
          character*1 ansrun
    c     Calculate the pseudorotation angle defined by atoms i1-inmem,
    c     Check if atoms ix(1-nmem) form a loop
              cospsrs(i)=0.d0
    c       As defined by Altona & Sundaralingam, JACS (1972) 94, 8205.
    c       Average over all cyclic permutations
    c        write (iout,8711) (theta(it)*radtodeg,it=1,5)
    c8711    format(' theta 1-5=',5f8.3)
    c     Cremer & Pople, JACS 97, 1354 (1975)
          call zeroit(r0,3)
          call zeroit(r1,3)
          call zeroit(r2,3)
          call vprd(r1,r2,rnorm)
    c      write (iout,2000) rnorm
    c2000  format(' rnorm=',3f10.6)
          cossum=0.0
              cossum=cossum+z(i)*cos(2*pi*m*(i-1-incgen)/float(nmem))
            cs=1.0/sqdenom
              cospsrs(m)=cospsrs(m)+cs
    c         write (iout,1004) psr5
    c1004  format(' Torsion-based pseudorotation angles with all origins:',
    c     -  /,5x,5f8.3,' deg')
          character* 132 line(maxrec)
          character*1 ansrun
          character*2 s5name(5),pname(5),s6name(6)
          character*4 atomnam
          character*15 question
            call getint('Number of atoms in the ring',27,nmem,1,MAXRING,
              call getint(question,15,999999,1,n,ix5(i),00)
            call getint('Residue number of the ring',26,0,1,numres,iresring,
    c       Find hexose sugar
    c       Apex is O
                call leftadjust4(atomnam,atomnam)
    c       Find pentose sugar
    c       Apex is O
                call leftadjust4(atomnam,atomnam)
    c       Find proline
    c       Apex is N (?)
    c           print *,'ia,r=',ia,line(index(ia))(inamcol1:inamcol1+1)
                call leftadjust4(atomnam,atomnam)
    c     Orient molecule along the diagonal
          closorgd=distminimg(cn,ih,n,nnh,edge,ioppbc,cell,ncell,ixyzhex,
    c     Calculate rotation matrix
    c     call chkort(rot)
          call rotate_c(cn,n,rot,cn,'DIAGONAL',8)
          closest=distminimg(cn,ih,n,nnh,edge,ioppbc,cell,ncell,ixyzhex,
            call extension(c,ih,nnh,1,nslt,cmin,cmax,c0,0,0,v)
              call genimdist(c(1,ia),cell,ncell1,ncell,icmin,rmin2)
          character*(*) label
    c     Find the orientation minimizing the size of the bounding cube
    c       3-D optimzation - initialize simplex
            call zeroit(simplex,12)
    c       print *,'-Initial simplex distances=',vertex
    c       call chkort(rot)
            call amoeba(simplex,vertex,4,3,3,ftol,iter,co,cnew,ih,nslt,nnh,
    c       2-D optimzation
            call opt2d(co,cnew,ih,n,nslt,nnh,edge,ioppbc,cell,ncell,
          call rotate_c(co,n,rot,cn,'OPTBOUND',8)
    c     Find the orientation maximizing image-image distances
    c     print *,'OPTIMIZE start closlim,n,nslt=',closlim,n,nslt
          call checkwall(ioppbc,edge,cmin,cmax,co,cell,1,ncell,
    c     print *,'edge=',edge,' ncell=',ncell
          closorg=distminimg(co,ih,nslt,nnh,edge,ioppbc,cell,ncell,ixyzhex,
            clo=distminck(co,ih,nslt,nnh,cell,cellalt,ioppbc,ncell,
    c       3-D optimzation - initialize simplex
            call zeroit(simplex,12)
    c       print *,'-Initial simplex distances=',vertex
    c       call chkort(rot)
            call amoeba(simplex,vertex,4,3,3,ftol,iter,co,cnew,ih,nslt,nnh,
    c       2-D optimzation
            call opt2d(co,cnew,ih,n,nslt,nnh,edge,ioppbc,cell,ncell,
          call rotate_c(co,n,rot,cn,'OPTIMIZE',8)
          closest=distminimg(cn,ih,nslt,nnh,edge,ioppbc,cell,ncell,ixyzhex,
              call cellreduce(0,cn,ih,nslt,nnh,edge,edgen,edge_gen,ioppbc,
              call trnsfr(edgenn,edgen,3)
              call trnsfr(edgenn,edge,3)
            call cellreduce(1,cn,ih,nslt,nnh,edgenn,edgen,edge_gen,ioppbc,
    c       Recreate the cell centers for edge
            call crorgn(edge,edge_gen,ioppbc,3,ncell,cell,cellalt,
    c     print *,'Reoptimize n,ioppbc=',n,ioppbc
          call optimize(c,c,c2,ih,nslt,nnh,n,edge,edgen,edge_gen,ioppbc,
    c       Select smallest volume
              closopt=closest
              call trnsfr(copt,c,3*n)
              call trnsfr(edgeopt,edgen,3)
              call matprod(rotopt,rot0,bestrot)
    c       Select largest minimum distance
              closopt=closest
              call trnsfr(copt,c,3*n)
              call matprod(rotopt,rot0,bestrot)
    c     Find the closest symmetrical image of the optimized ort to the original
            call unitmat(rot(1,1,k))
                call unitmat(rot0)
                    call matprod(rot(1,1,k),rot0,rot0)
                call rotate_c(c,n,rot0,c1,'SYMFIT1',7)
    c           call checkchir(c1,n,2,1,3,4,6,isg)
    c           print *,'k1,k2,k3=',k1,k2,k3,' isg=',isg
                  call trnsfr(rotopt,rot0,9)
    c              write (6,1781) ncha,k,d0,d1,rot0
    c1781          format(' nch,k=',2i3,' d0,d1=',2e12.5,' rot0=',(/,3f8.4))
            call rotate_c(c,n,rotopt,c,'SYMFIT',6)
            call matprod(rotopt,bestrot,bestrot)
          character*1 xyz(3)
          call zeroit(permute,9)
          call rotate_c(c,n,permute,c,'ALIGNLONG',9)
          call rotate_c(edge,1,permute,edge,'ALIGNLONG1',10)
          call matprod(permute,rot,rot)
          character*(*) label
          character*8 resnames(maxrsd)
          character* 132 line(maxrec)
          character*8 rnprev
          character*6 rnuprev
          character*19 segline
    c     Check if residue charges are integers, print charge sum
    c     print *,'CHECKRESCHARGESUM isegcol1,2=',isegcol1,isegcol2
    c       print *,'isg=',isg,' mlim=',molsltlim(1,isg),molsltlim(2,isg)
    c           resnames(iresprev)=rnprev
    c                   Divide extra charge (qdistr) into small chunks
                          call putreal(line(index(ja))(iqcol1:iqcol2),
                          charge(ja)=qfix
          character* 132 line(maxrec)
          common /logging/ logfile,ipredict
    c     Change atomnames using A for aromatic carbons to C
            call nextchar(line(index(ia)),ifc,inamcol2)
                  call askyn('Do you want to change A*** atoms to C***',40,
          character*8 resnames(numres),atnames(n),newname
          character*132 line(maxrec)
    c     print *,'CHECKFORDUP NDUPDEL,IDELDUP,IADDNUM,IASKDUP=',
    c    -  ndupdel,ideldup,iaddnum,iaskdup
                    call lastchar(atnames(ia),lc,nnamcol)
                      call askyn('Do you want to delete the duplicate',35,1,
    c                     Try adding numbers until nothing matches
                          call askyn(
                      call askyn(
                        call writeint(newname,lc+1,iadd,nnamcol-lc)
          character*1 ild(maxrsd)
          character* 132 line(maxrec)
          character*8 atnam
    c     print *,'CHPROTC n,maxrec,iverb,maxneig=',n,maxrec,iverb,maxneig
            call leftadjustn(atnam,atnam,8)
                call leftadjustn(atnam,atnam,8)
                call lookforneig(atnam(1:4),'C   ',ing,ic,nfound,ica,ierr)
                call lookforneig(atnam(1:4),'N   ',ing,in,nfound,ica,ierr)
                call lookforneig(atnam(1:4),'CB  ',ing,icb,nfound,ica,ierr)
                call lookforneig(atnam(1:4),'HA  ',ing,iha,nfound,ica,ierr)
    c             HA is missing only - generate it
                    c(k,maxrec)=c(k,ica)+
                  call checkchir(c,maxrec,ica,in,iha,icb,ic,isg)
    c           Non AA residue
    c           Both L and D was found
                call print1charlist(ild,ncheck,line,index,ifres,irescol1,
          character*4 atnam4,label
          character*1 charlist(nlist)
          character*(*) line(maxrec)
          character*1 aa1(50)
          character*8 resname
    c       Get 1-character AA list
              call leftadjustline(resname,1,lresname)
              call changeprot(resname,aa1(ii),2)
    c          if (i0 .eq. 0)
    c     -      write (6,7711) i,ii,index(i),lresname,resname,aa1(ii)
    c7711  format('  i,ii=',2i4,' index=',i4,' lresname=',i4,
    c     -  ' resname=',a,' aa1=',a)
    c     Check the chirality of atom c(1,i0). If i1,i2,i3,i4 are the indices of
    c     atoms with increasing atomic number then isg = +1 -> R ; isg = -1 -> S.
          call arrdiff(c(1,i1),c(1,i0),d01,3)
          call arrdiff(c(1,i2),c(1,i0),d02,3)
          call arrdiff(c(1,i4),c(1,i3),d34,3)
          call vprd(d01,d02,rnorm12)
    c     When it=0, find out if any of the cell sizes can be further reduced
    c     without affecting the closest approach
    c     When it=1, scale down the cell until the minimum approach is reached
    c     print *,'CELLRED  edge=',edge
    c     print *,'CELLRED  closmn=',closmn
            call crorgn(edge,edge_gen,ioppbc,3,ncell,cell,cellalt,
            closmin=distminimg(c,ih,n,nnh,edge,ioppbc,cell,ncell,ixyzhex,
            closmin=closmn
          call trnsfr(edgen,edge,3)
          call trnsfr(edgenn,edge,3)
          call zeroit(delxyz,3)
            call crorgn(edgenn,edge_gen,ioppbc,3,ncell,cell,cellalt,
    c       print *,'Cellred ic,closmin,dnew=',ic,closmin,dnew
    c         Successful reduction - see how far it can be taken
    c           write (6,8788) edgenn
    c8788        format(' edgenn=',3f10.5)
                call crorgn(edgenn,edge_gen,ioppbc,3,ncell,cell,cellalt,
    c             print *,'dnew,del=',dnew,del
                    call trnsfr(edgen,edgenn,3)
          call prtcell(ioppbc,edgen,edge_gen,0.0,volnew,nwnew,0)
          character*1 xyz
          common /axislab/ xyz(3)
    c     Calculate the cell parameters from the molecule size and the solvent
    c     layer width
    c       Mostly spherical PBC's (FCC, TO or HEX CP)
            call compact(co,cn,iatnum,ih,n,nnh,c0,rmin,rorgext,rorgcom,
    c         FCC
    c         Truncated octahedron
    c         edge parameter is the distance of the truncating squares
    c         from the center
    c         Inscribed sphere is tangent to centers of opposite hexagons
    c         Hexagonal close packing
    c       Rectangular, cubic or hexagonal
            call extension(co,ih,nnh,1,n,cmin,cmax,c0,0,1,v)
            call shiftmol(co,n,c0,cn,-1.0)
          call prtcell(ioppbc,edge,edge_gen,0.0,volnew,nwnew,1)
              call cellreduce(0,cn,ih,n,nnh,edge,edgen,edge_gen,ioppbc,npbc,
              call trnsfr(edge,edgen,3)
            call cellreduce(1,cn,ih,n,nnh,edge,edgen,edge_gen,ioppbc,npbc,
            call trnsfr(edge,edgen,3)
          common /numbers/ sq3,sq3inv,sq3p2,sq2p3
              cx_c=2.0*edge(1)
              cx_a=cx_c*sq3p2
            call vprd(edge_gen(1,2),edge_gen(1,3),c0)
    c#    MMC routine 064 lstmod: 09/12/90
    c*****Generate the negative cell-center coordinates for the pbc routines
          common /pbcrotmat/ torot_ac(3,3),torot_ca(3,3),tofac_ac,tofac_ca
    c     print *,'CRORGN ioppbc=',ioppbc,' edg=',edg
    c       Cubic
    c       Rectangular
    c       Face-centered cubic
    c       Hexagonal prizm (regular)
    c       Hexagonal prizm (skewed)
    c       Truncated octahedron
    cDAN    edge=edg(1)*2.0*sqrt(2.0)
    c       Hexagonal close packed
    c       Octahedral
            call zeroit(cx,3)
            cx(1)=edg(1)
    c       do k=1,3
    c         cy(k)=cx(k)*onethird-edge_gen(k,2)*sq22p3
    c       end do
    c       solve (ez.ex)=e**2/3; (ez.ey)=-e**2/3; (ez.ez)=e**2
    c       cz(3)=edg(1)
    c       cz(1)=(edg(1)**2/3.0*(cy(2)+cx(2))+edg(1)*(cy(3)-cx(3)))/
    c    -    (cx(1)*cy(2)-cy(1)*cx(2))
    c       cz(2)=(-edg(1)**2/3.0-cz(1)*cy(1)-cz(3)*cy(3))/cy(2)
    c       csum=cz(1)**2+cz(2)**2+cz(3)**2
    c       fac=edg(1)/sqrt(csum)
    c       do k=1,3
    c         cz(k)=cz(k)*fac
    c       end do
            cy(1)=-edg(1)/3.0
            cy(2)=edg(1)*sq22p3
            cy(3)=0.0
            cz(1)=edg(1)/3.0
            cz(2)=cy(2)/2.0
            cz(3)=edg(1)*sq2p3
            call trnsfr(edge_gen(1,1),cx,3)
            call trnsfr(edge_gen(1,2),cy,3)
            call trnsfr(edge_gen(1,3),cz,3)
    c        write (6,9877) edge_gen
    c       Sphere
    c     cell will contain the negative of the respective cell center coords
          call zeroit(cell,3*27)
    c     print *,'IOPPBC=',ioppbc
    c-----Siple cubic or rectangular
                  cell(j,ii)=x(j)
    c-----Face-centered cubic
            cell(2,i)=-edge
            cell(3,i)=-edge
          cell(2,2)=edge
          cell(3,2)=edge
          cell(2,4)=edge
          cell(3,5)=edge
            cell(1,4+i)=cell(2,I)
            cell(3,4+i)=cell(3,i)
            cell(1,8+i)=cell(2,I)
            cell(2,8+i)=cell(3,i)
          cell(1,14)=-2.0*edge
          cell(1,15)=-cell(1,14)
          cell(2,16)=cell(1,14)
          cell(2,17)=-cell(2,16)
          cell(3,18)=cell(1,14)
          cell(3,19)=-cell(3,18)
    c-----Hexagonal prism (regular)
    c     Vertex of the hexagon along the ixyzhex(2) axis
          cell(ixyzhex(3),2)=-edge*sq3p2
          cell(ixyzhex(2),2)=-edge*threp2
          cell(ixyzhex(3),3)=-edge*sq3
          cell(ixyzhex(3),4)=-edge*sq3p2
          cell(ixyzhex(2),4)=edge*threp2
            cell(ixyzhex(2),k)=cell(ixyzhex(2),l)
            cell(ixyzhex(3),k)=-cell(ixyzhex(3),l)
            cell(ixyzhex(1),l1)=edgex
            cell(ixyzhex(1),l2)=-edgex
            cell(ixyzhex(2),l1)=cell(ixyzhex(2),k)
            cell(ixyzhex(2),l2)=cell(ixyzhex(2),k)
            cell(ixyzhex(3),l1)=cell(ixyzhex(3),k)
            cell(ixyzhex(3),l2)=cell(ixyzhex(3),k)
    c-----Hexagonal prism (skewed)
    c     Vertex of the hexagon along the ixyzhex(2) axis
          cell(ixyzhex(3),2)=-edgex/2.0
          cell(ixyzhex(2),2)=-w
          cell(ixyzhex(3),3)=-edgex
          cell(ixyzhex(3),4)=cell(ixyzhex(3),2)
          cell(ixyzhex(2),4)=-cell(ixyzhex(2),2)
            cell(ixyzhex(2),k)=cell(ixyzhex(2),l)
            cell(ixyzhex(3),k)=-cell(ixyzhex(3),l)
            cell(ixyzhex(1),l1)=edgep
            cell(ixyzhex(1),l2)=-edgep
            cell(ixyzhex(2),l1)=cell(ixyzhex(2),k)
            cell(ixyzhex(2),l2)=cell(ixyzhex(2),k)
            cell(ixyzhex(3),l1)=cell(ixyzhex(3),k)
            cell(ixyzhex(3),l2)=cell(ixyzhex(3),k)
    c-----Truncated octahedon
    c     Truncated face transforms
    c     print *,'CRORG edge=',edge,' ncell=',ncell
            cell(i/2,i) =2.0*edge
            cell(i/2,i+1) = -2.0*edge
    c     +/-z, xy face transforms
                cell(1,ic)=(-1)**(ix)*edge
                cell(2,ic)=(-1)**(iy)*edge
                cell(3,ic)=(-1)**(iz)*edge
            call trnsfr(cellalt,cell,ncell*3)
            call rotate_c(cell,ncell,torot_ca,cell,'CRORGN',6)
    c-----Hexagonal close packed
    c     Neighbour cells in the hexagonal plane
          cell(1,2)=d
          cell(1,3)=d/2.0
          cell(2,3)=d*sq3p2
          cell(1,4)=d/2.0
          cell(2,4)=-d*sq3p2
          cell(1,5)=-d
          cell(1,6)=-d/2.0
          cell(2,6)=d*sq3p2
          cell(1,7)=-d/2.0
          cell(2,7)=-d*sq3p2
    c     Neighbour cells above with touching face
          cell(2,8)=d/sq3
          cell(1,9)=d/2.0
          cell(2,9)=-d/(2.0*sq3)
          cell(1,10)=-d/2.0
          cell(2,10)=-d/(2.0*sq3)
    c     Neighbour cells above with touching vertex
          cell(2,11)=-d*2.0/sq3
          cell(1,12)=d
          cell(2,12)=+d/sq3
          cell(1,13)=-d
          cell(2,13)=d/sq3
          call trnsfr(cell(1,14),cell(1,8),18)
            cell(1,k+6)=cell(1,k)
            cell(2,k+6)=cell(2,k)
            cell(3,k)=d*sq2p3
            cell(3,k+6)=-cell(3,k)
    c-----Octahedral
    c      write (6,9877) edge_gen
    c9877  format(' E_G: ',3(3f8.4,1x))
                    cell(k,ncell)=ix*edge_gen(k,1)+iy*edge_gen(k,2)+
            cz(k)=(edge_gen(k,1)+edge_gen(k,2)+edge_gen(k,3))/2.0
    c-----Input cell centers
    c-----Sphere
          common /boxdat/ xtlabc(6),xtlabc0(6),box(3),box0(3),edge_gen(3,3),
          character*1 xyz
          common /axislab/ xyz(3)
    c       For skewed hexagons, proportionality does not hold
              call trnsfi(ixyzhextraj,ixyzhextrajcharmm,3)
            call crorgn(edge,edge_gen,ioppbc,3,ncell,cell,cellalt,
                cellfac(k)=xtlabc(ixcrd(k))/xtlabc0(ixcrd(k))
                cellfac(k)=box(k)/box0(k)
                cell(k,ic)=cell0(k,ic)*cellfac(k)
    c      write (6,7671) xtlabc0,xtlabc,
    c     -  ((cell0(k,ii),k=1,3),ii=ncell-1,ncell),
    c     -  ((cell(k,ii),k=1,3),ii=ncell-1,ncell)
    c7671   format(' xtlabc0=',6f12.8,/,' xtlabc=',6f12.8,/,
    c     -  ' cell0i',6f12.8,' cell=',6f12.8)
    c     print *,'GENIMDIST ncell,ncell1=',ncell,ncell1
    c     print *,'GENIMDIST123 ncell,ncell1,ixyzexcld,ixyzincld=',
    c    -  ncell,ncell1,ixyzexcld,ixyzincld
            call genimdist(crm,cell,ncell1,ncell,icmin,rmin2)
    c       2D PBC
    c         Compare only with cells in the 2D plane
    c       1D PBC
    c         Compare only with cells along the 1D axis
    c     print *,'ICMIN=',icmin
    c      write (77,1717) ncell1,ncell,crm,(cell(k,1),k=1,3),rmin2
    c1717  format(' ncell1,ncell=',2i3,' crm=',3f10.5,' cell=',3f10.5,
    c    -  ' rmin2=',f10.5)
          common /distmindata/ ef1,ef2,ef3,ef4,edgep,edgep2,edgey,edgex,
          common /numbers/ sq3,sq3inv,sq3p2,sq2p3
    c       Face-centered cubic
    c       Regular hexagonal prizm
    c       Skewed hexagonal prizm
    c       Truncated octahedron
    c       Cubic
    c       Hexagonal close-packed
    c       Rectangular
          common /distmindata/ ef1,ef2,ef3,ef4,edgep,edgep2,edgey,edgex,
          common /numbers/ sq3,sq3inv,sq3p2,sq2p3
    c       Face-centered cubic
    c       Hexagonal prizm, prism along ixyzhex(1) ax, vertex on ixyzhex(2) ax
    c         Cell 6 or right side of cell 5
    c           Cell 5
    c           Cell 6
    c         Cell 1 or left side of cell 5
    c           Cell 5
    c       Rectangular
    c       Truncated octahedron, axes normal to the square faces
    c       Truncated octahedron, X axis normal to a hexagon face
            call genimdist(d,cell,1,ncell,icmin,rmin2)
    c       Hexagonal close-packed
    c       Distances from hexagonal plane neighbours
    c       Distances from face-touching upper neighbours
    c       Distances from vertex-touching upper neighbours
    c       Input pbc
            call genimdist(d,cell,1,ncell,icmin,sum)
    c       Check
    c       d(1)=dc1
    c       d(2)=dc2
    c       d(3)=dc3
    c       call genimdist(d,cell,1,ncell,icmin,sumx)
    c       if (abs(sum-sumx) .gt. 0.1 .or.
    c    -           icmin .eq. 1 .and. icell  .gt. 1 .or.
    c    -           icmin .gt. 1 .and. icell  .eq. 1)
    c    -     print *,'DMC ERROR: sum,x=',sum,sumx,
    c    -    ' icmin,icell=',icmin,icell
          character*1 xyz
          common /axislab/ xyz(3)
    c     Find an approximate center first
    c     write (06,*) 'EXTENSION nstart,n=',nstart,n,' nnh=',nnh
            cmin(k)=1.e+25
            cmax(k)=-1.e+25
    c           if (iprint .gt. 1) write (77,7711) ii,i,k,c(k,i)
    c           if (iprint .gt. 1) write (77,7711) i,i,k,c(k,i)
            c0(k)=(cmax(k)+cmin(k))/2.0
    c7711    format(2i5,' c',i1,'=',f10.5)
              call askstop(0)
    c     print *,'SHIFTMOL n=',n,' fac=',fac,' c0=',c0
              cnew(k,i)=c(k,i)+fac*c0(k)
          character*(*) lab
          common /rotwarn/ nrotwarn,nrottest,devmax
    c     Check the rotation matrix
    c           print *,'STOPPING'
    c           STOP
          character*(*) lab
          call check_rotmat(rot,lab,llab,ifail,iverb)
    c     Orient the molecule
            cn(1)=rot(1,1)*c(1,im)+rot(1,2)*c(2,im)+rot(1,3)*c(3,im)
            cn(2)=rot(2,1)*c(1,im)+rot(2,2)*c(2,im)+rot(2,3)*c(3,im)
            cn(3)=rot(3,1)*c(1,im)+rot(3,2)*c(2,im)+rot(3,3)*c(3,im)
            call trnsfr(cnew(1,im),cn,3)
    c       cnew(1,im)=cn1
    c       cnew(2,im)=cn2
    c       cnew(3,im)=cn3
          character*(*) q1,q2
          character*1 xyz(3)
          character*80 quiz
            call getreal(quiz,lquiz,default,c(k),noneg,ihelp)
          character*8 atnames(mxat),resnames(maxrsd)
          character*1 xyz,ans
          common /axislab/ xyz(3)
    c     print *,'TRANSROT maxrsd,mxat=',maxrsd,mxat
              call quiz(ans,itranstyp,' ',' ',0,
    c         Shift (first)
                call getxyz('Value to add to the ',20,' coordinates',12,0.0,
                call getint('Index of atom to be shifted',27,1,1,n,ias,0)
                call getxyz('New ',4,' coordinate of the selected atom',
                call quiz(ans,ictp,'g',' ',0,'molecular center',16,2,5,6,25)
                   call trnsfr(xyzmin,c0,3)
                  call cofms(c,xyzmin,n,atw)
              call shiftmol(c,n,xyzmin,c,shiftfac)
    c         Now reset into the PBC cell
              call setpbccell(' ',0,edge,edge_gen,cell,ncell,cellalt,
              call setrepats(nmolslt,molsltlim,nslt,nneig,ineig,
              call setpbcdim(ioppbc,ixyzhex,ixyzexcld,ixyzincld,xyz)
    c           Centering was requested
    c           Set atwtmp to zero for molecular ions
                call systemcenter(n,nmolslt,nmolsltnoion,molsltlim,c,ct,it1,
    c           Reset after shift
                    call cofms(c(1,molsltlim(1,is)),crep,natoms,atw)
                    call trnsfr(crep,c(1,molsltlim(3,is)),3)
                  call pbcreset(c(1,molsltlim(1,is)),natoms,crep,
                  call pbcreset(c(1,nslt+(iw-1)*nslv+1),nslv,
    c         Rotate by an input angle
              call genrot(rot,pi,iax,angle)
                call unitmat(rot)
                call rotate_c(c,n,rot,c,'TRANSROT',8)
              call getint('Atom index to put at the origin',31,1,1,nslt,
              call getint('Atom index to put on the X axis',31,1,1,nslt,
              call getint('Atom index to put in the X-Y plane',34,1,1,nslt,
              call trnsfr(xyzmin,c(1,iacent),3)
              call shiftmol(c,n,xyzmin,c,-1.0)
              call trnsfr(x,c(1,iaxax),3)
              call norm(x,1.0)
              call trnsfr(y,c(1,iaxyplane),3)
              call vprd(x,y,z)
              call vprd(z,x,y)
              call norm(y,1.0)
              call norm(z,1.0)
              call rotate_c(c,n,rot,c,'TRANSROTb',9)
              call getreal('Factor to multiply the coordinates with',39,
                  c(k,ia)=scfac*c(k,ia)
    c         Separate molecular complexes
                call askstop(1)
                  call getint('Index of solute molecule to move',32,1,1,
                  call getint('Index of solute molecule to move away from',
                call getreal('Distance to move',16,5.0,dmove,0,000)
                call setrepats(nmolslt,molsltlim,nslt,nneig,ineig,
                call zeroit(cent,3)
                      call trnsfr(ct(1,im),c(1,molsltlim(3,im)),3)
                      call extension(c,it1,0,molsltlim(1,im),
                      call zeroit(c00,3)
                          c00(k)=c00(k)+c(k,ia)
                        ct(k,im)=c00(k)/(molsltlim(2,im)-molsltlim(1,im)+1)
                      cent(k)=cent(k)+ct(k,im)
                  cent(k)=cent(k)/2.0
                    c(k,ia)=c(k,ia)+(y(k)-ct(k,im1))
    c         Swap chirality
              call getint('Index of the chiral center atom',31,999999,1,
              call getint('Index of the 1st neighbor to swap (1/2/3/4/)',
                call getint('Index of the 2nd neighbor to swap (1/2/3/4/)',
              call trnsfi(it1,nneig,nslt)
              call findbranch(iachir,iachir1,iachir3,c,nbranch1,ia_branch,
              call findbranch(iachir,iachir2,iachir3,c,nbranch2,
    c           Swap the coordinates of iachir1 and iachir2
                call unitvec(c(1,iachir),c(1,iachir1),x)
                call unitvec(c(1,iachir),c(1,iachir2),y)
                  c(k,iachir1)=c(k,iachir)+y(k)*rchir1
                  c(k,iachir2)=c(k,iachir)+x(k)*rchir2
                call buildbranch(c,nbranch1,ia_branch,
                call buildbranch(c,nbranch2,
          character*8 atnames(maxat),resnames(maxrsd)
    c     Atom ia_branch(i) has predecessors ia_predef(*,i),with distance,
    c     angle and torsion r_predef(i),ba_predef(i),ta_predef(i)
    c     move iachir to be the last neighbor of inchir
    c           write (6,2001) ianext,atnames(ianext)(1:latnam),
    c    -        resnames(ixres(ianext))(1:lresnam),iresno(ianext)
    c           write (6,2001) iatry,atnames(iatry)(1:latnam),
    c    -        resnames(ixres(iatry))(1:lresnam),iresno(iatry)
    c             write (6,2001) ip,atnames(ip)(1:latnam),
    c    -          resnames(ixres(ip))(1:lresnam),iresno(ip)
    c           write (6,2001) inchir,atnames(inchir)(1:latnam),
    c    -        resnames(ixres(inchir))(1:lresnam),iresno(inchir)
    c           Open ring, restart
                call getint('Bond # to break (1/2/3/...)',27,nmem-1,1,
                call trnsfi(nneig,nneig_sav,nslt)
                call breakbond0(i1,i2,nslt,nneig,ineig,ifail,maxneig)
                  call trnsfi(nneig,nneig_sav,nslt)
    c           print *,'Atom ',ianext,' added to ia_branch; nbr=',nbranch
    c     List is complete - set its properties
    c       Don't use the first atom on the list - it is swapped separately
    c       print *,'IAN,IA1,2,3=',ian,ia1,ia2,ia3
    c       print *,'R,BA,TA=',r_predef(iaa),ba_predef(iaa),ta_predef(iaa)
          call trnsfi(nneig,nneig_sav,nslt)
    c2001  format(' Ring member:',i6,1x,a,1x,a,i5)
    c       print *,'IA,IA1,2,3=',ia,ia_predef(1,ia),ia_predef(2,ia),
    c    -    ia_predef(3,ia)
            call addatom(1,c(1,ia_predef(1,ia)),c(1,ia_predef(3,ia)),
          character*1 xyz(3),ans
            call getint('Number of PBC dimensions',24,3,1,3,npbcdim,110)
                call quiz(ans,ixyzexcld,' ',
                call quiz(ans,ixyzincld,' ',
    c#    MMC routine 097 lstmod: 10/10/90
    c*****Generate an orientation matrix from the 3 euler angles
    c     this subroutine prepares a rotation matrix
    c     eq(4-47) of Goldstein (r=a)
          cfi=cos(fi)
          cps=cos(ps)
          cth=cos(th)
    c#    MMC routine 097 lstmod: 10/10/90
    c*****Generate an orientation matrix from the uniform distribution
    c     from the 3 Euler angles (i2dopt=0) or a random orientation, corresponding
    c     to a random rotation around axis i2dopt (i2dopt .gt. 0)
    c       eq(4-47) of Goldstein (r=a)
            call randpx(3,rn)
            cfi=cos(fir)
            cps=cos(psr)
            call randpx(1,rn)
            cfi=cos(fir)
            call zeroit(r,9)
          call arrdiff(r2,r1,e,3)
          call norm(e,1.0)
    c     Generate rotation matrix for 90 degree rotation around the ix axis
          call unitmat(r)
    c*****Set the seed of the random number generator (for reproducibility)
          common /rangen/ ixo
    c#    MMC routine 429 lstmod: 10/04/86
    c*****Congruential random number generator, Forsythe's constants
          common /rangen/ ixo
    c       Eliminate bits over 31
          call indexit(ix,1,n,0)
            call randpx(1,rand)
    c     Check rot for orthogonality
    c     Normalize d
    c*****c=a+b
            c(i)=a(i)+b(i)
    c*****c=a-b
            c(i)=a(i)-b(i)
    c*****c=a-b
            c(i)=a(i)-b(i)
    c#    MMC routine 277 lstmod: 03/29/00
    c*****dist2 = sum(ai-bi)**2
    c*****Real array transfer
    c*****Fast integer array transfer
    c*****Fast array transfer
    c*****Copy the bitpattern between an integer and a real
    c       Put the real into the integer
            call copyirr(realv,0)
            call copyiri(intgv,1)
    c       Put the integer into the real
            call copyiri(intgv,0)
            call copyirr(realv,1)
          common /janus/ realvar
          common /janus/ intgvar
          character*(*) line
    c     Finds the last nonblank in line
          character*1 char
          character*(*) line
    c     Finds the next character char in line
    c     print *,'FINDCHAR char=',char
          character*(*) line
    c     Finds the next nonblank in line
          character*1 tab,ctrlM
          common /tab/ tab,ctrlM
    c       if (line(i:i) .eq. ' ') print *,'nextc i=',i,' blank'
    c       if (line(i:i) .eq. tab) print *,'nextc i=',i,' tab'
    c         if (line(i:i) .eq. '!') ifc=len
          character*(*) line
    c     Finds the next blank in line
          character*1 tab,ctrlM
          common /tab/ tab,ctrlM
    c         if (line(i:i) .eq. '!') ifc=len
          character*(*) line
          call nextchar(line,ifc,len)
    c       Look for closing quote
              call findnextchar('"',line,ifc,len)
              call findnextchar("'",line,ifc,len)
            call nextblank(line,ifc,len)
          character*(*) line
            call nextstring(line,ic,icf,icl,len)
    c     print *,'--- 2D optimization started'
    c     print *,'Initial touch=',dl,' i2dopt=',i2dopt
    c     print *,'al,am,ar=',al,am,ar
    c     print *,'dl,dm,dr=',dl,dm,dr
    c         Middle is the lowest
    c         Left is the lowest
    c         Right is the lowest
    c     The molecule co will be replaced by a molecule having the same
    c     internal geometry as rc with the com's coinciding and the best
    c     orientational overlap
    c     Obtain  orientation of co
    c     print *,'fixup n=',n
    c     write (6,1002) co
    c1002  format(' co=',/,(3f10.5))
          call cofms(co,como,n,aw)
          call trnsfr(a,rc,9)
          call ormat(orient,a,b,n,1,6)
          call rotate_c(rc,n,orient,cn,'FIXUP',5)
          call shiftmol(cn,n,como,cn,1.0)
    c     Create a standard water
            c(i,1)=0.0
            c(3,i)=0.0
          c(1,2)=xw
          c(1,3)=xw
          c(2,2)=yw
          c(2,3)=-yw
    c#    MMC routine 017 lstmod: before 07/15/85
    c*****compute the c.o.m. centered coordinates from cs into csr
          call cofms(cs,cmcrd,natms,aw)
              csr(k,i)=cs(k,i)-cmcrd(k)
    c#    MMC routine 021 lstmod: 06/18/93
    c*****compute center of mass from global atomic coordinates
    c     print *,'COFMS natm=',natm
    c     print *,'COFMS aw=',aw
    c     write (6,1003) csa
    c1003  format(' csa=',/,(3f10.5))
    c#    MMC routine 022 lstmod: 12/08/86
    c     3-atom part originally written by P.K. Mehrotra.
    c      b=orm.a
    c     ++++++++++++++++++++++++++++++++++++++
    c     a - coordinates of the molecule before rotation ( local system)
    c     b - coordinates of the molecule after rotation (global syetm).
    c     It is assumed that the coordinates of the molecule are given
    c     columnwise, i.e., the coordinates of the atom 1 constitute
    c     the first column, coordinates of the atom 2 constitute the
    c     the second column, etc.
    c     Choose the first atom as the origin
    c     print *,'ORMAT iverb,natm=',iverb,natoms
    c     Generate three specific orthogonal vectors in lab frame
          call vprod(d,2,3,1)
          call vprod(d,1,2,3)
          call vprod(e,2,3,1)
          call vprod(e,1,2,3)
    c     Check for colinearity of the atoms
    c     switch to two-atom algorithm
    c     Diatomic
    c     If the two atoms coincide,generate unit matrix
    c         Zero component found
    c     No zero component
    c         Zero component found
    c     No zero component
          call mnorm(e)
    c     Now, d=orm*e, thus orm=d*inv(e) and
    c     The inverse of an orthonormal matrix is its transpose
          call trnsfr(r3,r,9)
          call trnsfrd(dout,d,3)
    c#    MMC routine 025 lstmod: 02/06/86
    c*****Computes the vector product of the columns i and j into k
    c#    MMC routine  61 lstmod: 02/06/86
    c*****Computes the vector product a x b and saves it into c
          c(1)=a(2)*b(3)-b(2)*a(3)
          c(2)=a(3)*b(1)-b(3)*a(1)
          c(3)=a(1)*b(2)-b(1)*a(2)
    c#    MMC routine 048  lstmod: 11/13/90
    c*****Computes the vector product a x b and saves it into c
          c(1)=a(2)*b(3)-b(2)*a(3)
          c(2)=a(3)*b(1)-b(3)*a(1)
          c(3)=a(1)*b(2)-b(1)*a(2)
    c#    MMC routine 052 lstmod: 11/13/90
    c*****Normalizes the matrix r
          call zeroiti(indexrev,0,max)
    c*****Calculate the angle c2-c1-c3 and distance c1-c2; c1-c3
          cosa=(rroh1+rroh2-rrhh)/(2.d0*droh1*droh2)
    c     if iaaconv=1 convert from 1-digit to 3 digit AA code
    c     if iaaconv=2 convert from 3-digit to 1 digit AA code
          character*8 resnam
          character*1 aanames1,resnam1
          character*2 mmodtoamb
          character*3 aanames3
          common /atnamcon/ mmodtoamb(100),aanames1(58),aanames3(58),
            call leftadjustn(resnam,resnam,8)
    c     Determine atomic number from an atomname
          character*8 atomnam
          character*1 anam1(20)
          character* 132 line
          call nextchar(line,ic,132)
              call zeroiti(ineig(1,i),0,maxng)
          character* 132 line(maxrec)
          character*4 atnami
          character*8 resnami
    c     Set up neighbour list
          character*4 namfcg
          common /connatdat/ ramax(99),ramax2(99),hlimfac,ianfg(99),
          common /clonedat/ nclone,iaclnf(1000),iaclnl(1000),ncopcln(1000)
    c     print *,'NNLIST n,nslt,islvw,ihbondcalc=',n,nslt,islvw,ihbondcalc
          call nninit(nneig,nhbneig,ineig,nhneig,nnneig,ncneig,
    c     Generate connectivity for solvents
    c     print *,'NNLIST n,nslt,nslv,numsolv=',n,nslt,nslv,numsolv
            call nnlist0o(nslt+(is-1)*nslv+1,nslt+is*nslv,iatnum,c,nneig,
            call nnlist0(1,nhb,nslt,islvw,iatnum,ifchrg,c,nneig,
            call nnlist0(1,iaclnf(1)-1,nslt,islvw,iatnum,ifchrg,c,nneig,
    c         Create list separately for the clones
              call nnlist0(ifirst,ilast,nslt,islvw,iatnum,ifchrg,c,nneig,
    c        Copy cloned nn info
                 call trnsfi(ineig(1,ianew),ineig(1,ia),maxng)
    c       Repeat for test
            call nninit(nneig,nhbneig,ineig,nhneig,nnneig,ncneig,
            call nnlist0o(1,ntest,iatnum,c,nneig,nhbneig,
    c     Save the pseudo-atom neighbours for the real atoms
    c     Functional group search will not see the pseudo-atom neighbours
    c         Pseudo atom found
    c             Add i to the neighbour list of ia
    c     Check for unconnected atoms not already made a molecule
          common /connatdat/ ramax(99),ramax2(99),hlimfac,ianfg(99),
          character* 132 line(maxrec)
    c     print *,'NNL0 n,nslt,islvw,ihbondcalc=',n,nslt,islvw,ihbondcalc
    c       LES structure, call nnlist00 with the whole range
            call nnlist00(nfirst,n,nslt,islvw,iatnum,ifchrg,c,nneig,nhbneig,
    c       Call the actual near neighbour search by segments
    c         Find limits of the next segment
    c           If segments are not contiguous, just call nnlist00 for
    c           the whole range
                  call nninit(nneig,nhbneig,ineig,nhneig,nnneig,ncneig,
                  call nnlist00(nfirst,n,nslt,islvw,iatnum,ifchrg,c,nneig,
              call nnlist00(ifat,ilat,nslt,islvw,iatnum,ifchrg,c,nneig,
          call nnlist0o(nfirst,n,iatnum,c,nneig,nhbneig,ineig,nhneig,nnneig,
          common /connatdat/ ramax(99),ramax2(99),hlimfac,ianfg(99),
          character* 132 line(maxrec)
          character*4 atnami,atnamj
          character*8 resnami,resnamj
    c     Set up neighbour list using linked cells
    c     nhneig, ncneig,nnneig,nsneig,npneig:
    c     Number of H, C, N, S, P neighbours, resp.
    c     print *,'NNLIST00 nfirst,n,nslt,islvw=',nfirst,n,nslt,islvw
    c     Find minimum and maximum atomic radius and extents of the molecule
    c       write (77,9877) i,line(index(i))(inamcol1:inamcol2),
    c    -    iatnum(i),ramax2(iatnum(i)),ra2max
    c9877   format(' NNLIST0',i6,a,' ia',i3,' ramax2=',f5.2,' ra2max=',f5.2)
    c     print *,'DIV=',div,' RAMX=',ramx
          call extension(c,nneig,0,nfirst,n,xyzmin,xyzmax,centinp,0,0,v)
          call gridspace(c,isegno,nfirst,n,xyzmin,xyzmax,div,
    c     Loop on the boxes
    c                   ni, nj are the number of atoms in the box (i1,i2,i3) and
    c                   (j1,j2,j3)
    c                       write (78,*) 'NNLIST00 i,j=',i,j,' numhyd=',numhyd
    c                         At most one hydrogen
    c                           Bond found. Bonds to electron, charge or lone pair
    c                           are not saved in the 'real' atoms' list here.
    c                             Solvent-solvent connectivity is already done
                                  call savebond(i,j,iatnum,nneig,ineig,
                                  call savebond(j,i,iatnum,nneig,ineig,
    c                           No bond, but may be hydrogen bond
    c                             No H-C or H-+
    c                             slv-slv bond only allowed for water bridge calc
    c??                           Hbond with C between slt-slv?? - removed!!
    c                               if (i .eq.  900 .and. j .eq. 3142 .or.
    c    -                              j .eq.  900 .and. i .eq. 3142)    
    c    -write (6,8976) ,i,j,r2,rlimhb2,ramax2(iatnum(i)),ramax2(iatnum(j))
    c8976 format('i,j=',2i6,' R2=',f6.2,'RLIMHB2=',f6.2,' ramax2:',2f7.2) 
                                  call maybehbond(r2,i,j,nneig,nhbneig,
    c                        print *,'MAYBE H-+:',i,j,nhbneig(i),nhbneig(j)
    c                        print *,'MAYBE H-+: r2=',r2,' rlimhb2=',rlimhb2
    c     Now screen the H-bonds for the angle and for C-H...X bond
    c         Atom ia is always the H of the H bond
    c         ihb0 is the heavy atom of the donor H
              call get_heavyat(ia,nneig,ineig,ixres,nframe,ihb0,maxng,
    c             ihb0 is the heavy atom of the donor H
    c             ihb0=0
    c             nng=nneig(ia)
    c             do while (nng .gt. 0)
    c               ihb0=ineig(nng,ia)
    c               if (iatnum(ineig(nng,ia)) .eq. 1)  then
    c                 ihb0=ineig(nng,ia)
    c                 nng=0
    c               end if
    c               nng=nng-1
    c             end do
                    call checkhbclose(c,n,ia,ihb0,
                    call checkhbclose(c,n,ihb0,ihb,
                    call angdistw(c(1,ia),c(1,ihb),c(1,ihb0),rHB,rb,rab,ang)
    c               Hydrogen had no heavy atom bonded to it - discard
    c               Remove ihb from the HB list of ia
    c               write (77,*) 'rHB,rb,rab=',rHB,rb,rab
    c               Remove ia from the HB list of ihb
    c         Cation - only H-bonds to water oxygen
    c             Not dropped yet
    c               write (77,*) 'ia,ihb,nh12=',ia,ihb,nh12
                    call angdistw(c(1,ihb),c(1,ia),c(1,nh12(1)),rHB,rb,
                      call angdistw(c(1,ihb),c(1,ia),c(1,nh12(2)),rHB,rb,
    c               Remove ia from the HB list of ihb
    c     Now condense the list to fill in the zeros
          call checkhblist(n,ineig,nhbneig,maxng)
    c     print *,'GET_HEAVYAT ia=',ia,' MAXNEIG,MAXREC=',maxneig,maxrec
          character*(*) bondlab
    c     Set up hydrophobic bond list using linked cells
    c     print *,'NNLISTHPH_SLTB n,nosameseg,iselfanc,rhphmax=',
    c    -  n,nosameseg,iselfanc,rhphmax
    c     indexa(ia)=1 or 2: Anchor atom
    c     indexa(ia)=-1 or -2: Non-anchor atom, but can form HPH bond/salt bridge
    c     iabs(indexa(ia))=1: +; =2: -
          call extension(c,nneig,0,1,n,xyzmin,xyzmax,centinp,0,1,v)
          call gridspace(c,isegno,1,n,xyzmin,xyzmax,div,
    c     All hydrophobic carbons are indexed in the grid
    c     Loop on the boxes
          call zeroiti(nneig,0,n)
          call zeroiti(nhbneig,0,n)
    c                   ni, nj are the number of atoms in the box (i1,i2,i3) and
    c                   (j1,j2,j3)
    c                         Salt bridge - exclude carbons
    c                         Make sure ia and ja far enough apart in topology
    c                         write (77,7211) ia,ja,idoit,
    c    -                      (ineig(in,ia),in=1,n14neig(ia))
    c7211                      format(' ia,ja,idoit=',3i5,' in(ia)=',6i5)
    c                           Bond found
    c     Set up mutually proximal pair list using linked cells
    c     print *,'NNLISTMPX n=',n,' nanchorr,nanchorn=',nanchorr,nanchorn
    c     print *,'NNLISTMPX n=',n,' maxbox,maxrec=',maxbox,maxrec
    c     indexa(ia): Reference set atoms
    c     indexov(ia): neighbour set atoms
          call extension(c,nneig,0,1,n,xyzmin,xyzmax,centinp,0,1,v)
          call gridspace(c,isegno,1,n,xyzmin,xyzmax,div,
          call zeroiti(it1,0,n)
    c     Leave hydrogens out from the mpx search
    c     Sort the indices into contiguous list NONE/ANCHOR/NEIGHBOUR (0/1/2)
            call indexit(it2,1,nb,0)
            call mrgsrt(6,it2,temp1,nb,it3,it4,it5,temp2,n)
            call trnsfi(it5,indices(1,ib),nb)
    c       write (88,6941) 'it5',ib,(it5(ia),ia=1,nb)
    c       write (88,6941) 'it1',ib,(it1(it5(ia)),ia=1,nb)
    c       write (88,6941) 'ian',ib,(iatnum(it5(ia)),ia=1,nb)
    c       write (88,6941) 'it2',ib,(it2(ia),ia=1,nb)
    c       write (88,6941) 'ind',ib,(indices(ia,ib),ia=1,nb)
    c6941   format(1x,a,i6,':',(20i6))
    c     Find the limits of the NONE and ANHOR stretches
    c             No ones
    c           All zeros
    c         write (88,6791)ib,it3(ib),it4(ib),
    c    -      (it1(indices(ia,ib)),ia=1,nbox(ib))
    c         write (88,6792)ib,(indices(ia,ib),iatnum(indices(ia,ib)),
    c    -      ia=1,nbox(ib))
    c6791     format(i6,' it3,it4=',2i4,' it1(indices)=',(10i6))
    c6792     format(i6,' indices,iatnos=',(10i6))
    c     Anchor atoms: ia=it3(ib),it4(ib)-1
    c     Neighbour atoms: ia=it4(ib),nbox(ib)
    c     All heavy atoms are indexed in the grid
    c     Loop on the boxes
          call zeroiti(nneig,0,n)
          call zeroiti(nhbneig,0,n)
          call zeroiti(it2,0,n)
    c     it2(ia), temp2(ia): nearest neighbor atom and distance from anchor ia
    c     and vice versa
    c                   ni, nj are the number of atoms in the box (i1,i2,i3) and
    c                   (j1,j2,j3)
    c                     Work with atom pairs in boxes indexi and indexj
    c9783     format(' ia=',i6,' it1,it2=',2i6,' r21,r22=',2f8.2)
    c     print *,'GRIDSPACE MAXBOX,MAXREC=',maxbox,maxrec,' DIV=',div
    c       Increase gridsize to reduce the number of boxes under maxrec
          call zeroiti(nbox,0,ngrid)
    c     print *,'NNLIST00 NFIRST,N=',nfirst,n
    c     nboxmax=0
    c         Save i into the box represented by the indices ix(1-3)
    c         if (nboxmax .lt. nbox(indexi)) nboxmax=nbox(indexi)
    c     print *,'NBOXMAX=',nboxmax
          character*4 atnami
          character*8 resnami
    c     write (78,*) 'MAYBEHBOND i,j=',i,j
    c         Hydrogen bond found
          character*(*) label
    c     Set up neighbour list using linked cells
    c     print *,'NNLS maxng,nnlistlen,maxbox=',maxng,nnlistlen,maxbox
          call extension(c,nneig,0,nfirst,n,xyzmin,xyzmax,centinp,0,0,v)
          call zeroiti(nneig,0,n)
    c       Increase gridsize to reduce the number of boxes under maxrec
          call zeroiti(nbox,0,ngrid)
    c        Save i into the box represented by the indices ix(1-3)
    c     print *,'nboxmax=',nboxmax
    c     Loop on the boxes
    c                   ni, nj are the number of atoms in the box (i1,i2,i3) and
    c                   (j1,j2,j3)
    c                         Bond found.
    c     print *,'maxnn=',maxnn
          common /connatdat/ ramax(99),ramax2(99),hlimfac,ianfg(99),
          character*4 atnamj
          character*8 resnamj
    c     Save atom i as the neighbor of atom j
    c           Keep the neighbor list sorted
          call zeroiti(nneig,nfirst-1,n)
              call decidebondcut(iatnum(i),iatnum(j),rlim)
          common /connatdat/ ramax(99),ramax2(99),hlimfac,ianfg(99),
          character* 132 line(maxrec)
          character*4 atnami
          character*8 resnami
    c     Set up neighbour list using the trivial algorithm
    c     nhneig, ncneig,nnneig,nsneig,npneig:
    c     Number of H, C, N, S, P neighbours, resp.
    c       write (77,7711) i,iatnum(i),ramax(iatnum(i)),ramax2(iatnum(i))
    c7711    format(i5,' iano=',i2,' ramax,2=',2f10.5)
              call decidebondcut(iatnum(i),iatnum(j),rlim)
    c           Bond found. Bonds to electron, charge or lone pair are not
    c           saved in the 'real' atoms' list here.
                call savebond(j,i,iatnum,nneig,ineig,nhbneig,nhneig,ncneig,
                call savebond(i,j,iatnum,nneig,ineig,nhbneig,nhneig,ncneig,
    c           No bond, check if hydrogen bond
    c               Hydrogen bond found
          character* 132 line(maxrec)
          character*2 iatnm2
          common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
    c     Check for atoms exceeding their valence
                call decidebondcut(iatnum(i),iatnum(j),rlim(jj))
    c     Check for bonded atoms far apart
          character* 132 line(maxrec)
          common /connatdat/ ramax(99),ramax2(99),hlimfac,ianfg(99),
          character*2 iatnm2
          common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
          character*4 pbclab
          character*4 atnami
          character*8 resnami
    c     Check solute topology
    c       Check bonds if they are too short
                  call decidebondcut(iatnum(i),iatnum(j),rlim)
    c               Too short bond found.
    c       Check for close approach
    c           Check solute intramolecular NB distances
                      call decidebondcut(iatnum(i),iatnum(j),rlim)
    c                   Contact between atoms found.
    c           Check solute intermolecular distances
    c                   Check for images
                        call arrdiff(c(1,i),c(1,j),rpbc,3)
                        call distmincalc(ioppbc,cell,ncell,ixyzhex,edge,
                      call decidebondcut(iatnum(i),iatnum(j),rlim)
    c                   Contact between atoms found.
                          call arrdiff(c(1,i),c(1,j),rij,3)
                          call genimdist(rij,cell,1,ncell,icminn,rmin2)
    c     Check for solute-solvent contacts
              call arrdiff(c(1,ia),c(1,jm0+1),rpbc,3)
    c           Check for images
                call distmincalc(ioppbc,cell,ncell,ixyzhex,edge,
    c           General image search to get neighbor cell no (could be improved)
                  call arrdiff(c(1,ia),c(1,jm0+ja),rpbc,3)
                call decidebondcut(iatnum(ia),iatnum(ja),rlim)
                    call arrdiff(c(1,ia),c(1,ja),rij,3)
                    call genimdist(rij,cell,1,ncell,icminn,rmin2)
    c     Check for solvent-solvent contacts
              call arrdiff(c(1,im0+1),c(1,jm0+1),rpbc,3)
    c           Check for images
                call distmincalc(ioppbc,cell,ncell,ixyzhex,edge,
                    call arrdiff(c(1,im0+ia),c(1,jm0+ja),rpbc,3)
                  call decidebondcut(iatnum(iac),iatnum(jac),rlim)
                      call arrdiff(c(1,iac),c(1,jac),rij,3)
                      call genimdist(rij,cell,1,ncell,icminn,rmin2)
          character*(*) label,bondtype,inpfile
          character*80 bond,listfile
          character*132 line(maxrec)
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
          common /bondpairs/ ihbpair(2,MAXBONDS),ihb_pair_res(3,MAXBONDS)
          character*6 bondlab
    c     print *,'FILTERBOND NBFOUND=',nbfound
          call getfiltlims(0.0,100.0,1,nresslt-1,percmin,percmax,
          call zeroit(percbond,n)
                  call bonddescr(it1(ib),ihbpair,line,index,iresno,isegno,
                  call bonddescr(it2(ib),ihbpair,line,index,iresno,isegno,
                  call bonddescr(it3(ib),ihbpair,line,index,iresno,isegno,
                  call bonddescr(it4(ib),ihbpair,line,index,iresno,isegno,
          call zeroiti(irrix,0,n)
          call zeroiti(itemp1,0,numres)
          call askyn('Do you want to calculate bond autocorrelation',45,
          call askyn(
            call changeext(inpfile,listfile,namleni,llistfile,'pls',3,0,0)
            call openfile(50,0,' ',1,'new',listfile,llistfile,notfnd,0,1,1,
            call bonddescr(i,ihbpair,line,index,iresno,isegno,inamcol1,
    c        write (06,9781) ib,i,nhbpers(i),nhbdist(i)
    c9781    format(' i=',i5,' ib=',i5,'  nhbpers=',i5,' nhbdist=',i5)
              call getbondtrack(i,itrack,ifirstframe,lastframe,30,nframe)
              call autocorr(ib,i,itrack,ifirstframe,lastframe,iframeunit,
          call checkdim(mxcopy+nauc_extra,mxcopy,'MAXCOPY',7,
          character*(*) label,bondtyp
    c     print *,'GETFILTLIMS label=',label,' bondtyp=',bondtyp
          call getreal('MINimum percentage presence',27,percmind,percmin,1,
          call getreal('MAXimum percentage presence',27,percmaxd,percmax,1,
          call getint('MINimum residue-residue sequence distance',41,
          call getint('MAXimum residue-residue sequence distance',41,
    c     Establish auc caculation details
          character*(*) label
          character*6 frunit
          common /frameunit/ lfrunit(4),frunit(4)
          character*1 ans
          character*60 auc_type(5)
          call quiz(ans,iauctype,'o',' ',0,
            call getint('Frame number to pad the tracks to',33,mxframes,1,
            call getint(
              call getint('Maximum number of frames to reuse',33,nframe,
              call getint(
          call askyn(
          character*132 line(maxrec)
          character*(*) bondlab
          character*2 sc_lab(2,2)
          character*(*) bond
          character*6 frunit
          common /frameunit/ lfrunit(4),frunit(4)
    c     write (iout,*) 'PERSISTENCE ITF,ITL=',itf,itl,
    c    -  ' IUSELASTON=',iuselaston
    c#    MMC routine  98/b lstmod: 06/21/05
    c*****List bonds, angles, torsions as requested
          character* 132 line(maxrec)
          character*8 lab
                cv(in)=sqrt(dist2(cslt(1,ia),cslt(1,ineig(in,ia))))
    c           Generate angle list
    c           Generate torsion angle list
    c                       ang(nang)=dihangl(cslt(1,inn1),cslt(1,ia),
    c    -                    cslt(1,inj),cslt(1,inn2),1,iout)
    c#    MMC routine 164 lstmod: 07/11/96
    c*****Assigns the appropriate fg types to the atoms in c
          character* 132 line(maxrec)
          character*2 iatnm2
          common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
          character*1 sp
          character*4 namfcg
          common /connatdat/ ramax(99),ramax2(99),hlimfac,ianfg(99),
    c     c(1,i) : coordinates of atom i
    c     n : number of atoms in c
    c     icl(i) : type af atom i
    c     ian(i) : atomic number of atom i
    c     nfg: number of functional groups found
    c     nfgmem(i): number of atoms in fcg i
    c     ifgstr(i): start of member list for fcg i in indxat
    c     itypfg(i): functional group type of atom i
    c     ifgtyp(i): the functional group type if the i-th fcg before sorting
    c     indxfg(i): the i-th atom belongs to the indxfg(i)-th fcg
    c     ifgaix(i): the original index of the i-th atom in the atomlist
    c     ixfg(ig) : the ig-th functional group after sorting
    c     sorted by functional groups
    c     nneig(i) : number of neighbours of atom i
    c     ineig(1,i) : list of neighbours of atom i
    c     nneigh(i) : number of hydrogen neighbours of atom i
          call zeroiti(itypfg,n0-1,n)
    c     Check for "valence errors"
            call askyn('Do you want to break any bond',29,1,1,ibreak,0,0)
                call getintline(
                call breakbond(in12(1),in12(2),n0,n,nneig,ineig,nneiga,
                call breakbond(in12(2),in12(1),n0,n,nneig,ineig,nneiga,
    c     if (nverr .gt. 0) go to 120
    c     Search for O and P first
    c         Oxygen found
    c          print *,'i,nneiga(i),ineig(1,i)=',i,nneiga(i),ineig(1,i)
    c             COO- found typ=27
    c             >C=O found,  typ=17
    c           Ester oxygen found, typ=19 or 20 (for phospho ester)
    c           -OH found, typ=21
    c         Phosphorus found  (>PO2 : ityp 22)
    c     Search for carbon and nitrogen next
    c         Unassigned carbon found
    c         Unassigned nitrogen found
    c         Four neighbours, assumed to be positively charged
    c     Label hydrogens on carbonyl
    c         Hydrogen on a C=O is type 18
    c     Search for -S- (type 28) and -SH (type 29)
    c         Sulphur found
    c           -S-, -SH or HSH found
    c             -SH or HSH found
    c             -S- found
    c     Assign the single atom funcional groups
    c           If no neighbours, must be an ion
    c     Label all unassigned atoms
    c     Sort atoms by functional groups
    c     index (sort) fcg's by type
    c      if (nfg .gt. 1) then
    cc       Sort functional groups
    c        do 50 i=1,nfg
    c          j1=i+1
    c          do 50 j=j1,nfg
    c            if (ifgtyp(i) .gt. ifgtyp(j)) then
    c              ii=ifgtyp(i)
    c              ifgtyp(i)=ifgtyp(j)
    c              ifgtyp(j)=ii
    c              ii=ifgstr(i)
    c              ifgstr(i)=ifgstr(j)
    c              ifgstr(j)=ii
    c              ii=nfgmem(i)
    c              nfgmem(i)=nfgmem(j)
    c              nfgmem(j)=ii
    c            end if
    c50      continue
    c      end if
    c       Print list
    c     print *,'BREAKBIND0 nneig(i1,i2),n=',nneig(i1),nneig(i2),n
    c     Establish the bond orders for Macromodel format (as best as possible)
          character* 132 line(maxrec)
          character*2 iatnm2
          common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
          character*4 dbonds(3,100),dbres(100),thisres,thisresl,
          character*8 resnam
    c     For residue dbres(1,i) there is a double bond between atoms
    c     dbres(2,i) and dbres(3,i)
    c     All names leftadjusted
    c     Scan dbond list to find a list of special residues
    c         New residue in the double bond list
    c     Find hybridization
    c        write (77,9857) i,iatnum(i),(ineig(j,i),j=1,nneig(i))
    c9857    format(i5,' atnum=',i3,' ing=',20i5)
    c             Calculate angle betwen in-th and (in+1)th neighbor
                  cosa=dble(d12s)/sqrt(d1s*d2s)
    c            write (199,7171) i,line(index(i))(inc1:inc2),angav,mmtype(i)
    c7171        format(i5,2x,a4,' angle=',f8.3,' mmtype(i)=',i2)
    c     See end-node carbons and nitrogens
    c     Find bond orders
    c       nvalx is the excess valence
    c         For carbons than might be united atoms, deduce united atoms
    c         Do end node atoms last
    c     Look at end-node carbons
    c     Initialize bond orders to one
    c     Scan the atom list to add the explicitly given double bonds
    c         New residue
              call leftadjust4(thisres,thisresl)
    c         See if there is double bond in the list for this residue
              call leftadjust4(line(index(i))(inc1:inc2),atnaml)
    c           Check if atom (i) is in the dblist
    c             Check the bonds of atom i against the list
                    call leftadjust4(line(index(in1))(inc1:inc2),atnam1l)
    c                 Match found
    c                 print *,'Added i,in1=',i,in1
    c     Take care of atoms with one neighbours first, in an iterative fashion
    c           Multiple terminal bond found, add to bond orders
    c           Drop atom from further consideration
                call swapng(ineig,ibnd,in1,jn1,nn1,nattot,maxneig)
    c             Immediate neighbour became an end node - repeat
    c     At this point only atoms in loops are left
    c         Distribute additional bond orders
    c         Find the neighbor with the lowest bond order
    c     Now see atoms with nvalmax > nval (like N)
    c             Check if any neighbor has unsatisfied valence
    c                 Increase bond order of i-i1)
    c       See if loop double bonds can be swapped around
    c          write (77,9571) i,(ineig(j,i),j=1,nneig(i))
    c9571      format(i4,' ing=',20i5)
    c           Scan neighbor for atoms with multiple bond
    c             write (77,*) 'i,in,i1=',i,in,i1
    c                write (77,*) 'i,in1,i2=',i,in1,i2
    c                 Multiple bond was found in the neighborhood: i1-i2
    c                 Check if i2 has a neighbor with unsatisfied valence
    c                   write (77,*) 'i,in2,i3=',i,in2,i3
    c                     Decrement bond order of i1-i2; increment i-i1 and i2-i3
    c                     First move the neig info to their new place
                          call swapng(ineig,ibnd,i2,in2,nntemp(i2),nattot,
    c                       write (77,*) 'i,i3,nntemp(i3)=',i,i3,nntemp(i3)
                          call swapng(ineig,ibnd,i,in1,nntemp(i),nattot,
    c      do i=1,nattot
    c        write (100,7711) i,iatnum(i),nneig(i),nvalx(i),
    c     -      (ineig(j,i),ibnd(j,i),j=1,nneig(i))
    c7711    format(i4,' atno=',i3,' nn=',i2,' nvalx=',i3,
    c     -    ' ineig,iibnd=',6(i4,i2,2x))
    c      end do
    c     Finally, set atomtypes
    c           Unsaturated carbon --> united atom
    c         Saturated carbon's atomtype is just the hybridization number
    c1201  format(1x,i3,6(1x,i5,1x,i1))
    c#    Based on MMC findcent
    c*****Finds the backbone of a molecule
          character* 132 line(maxrec),blankline,linep
    c     First grow neighbour list from n0
    c     print *,'FINDBACKBONE n0,n,maxbox=',n0,n,maxbox
            call growchain(n0,n,0,ifirst,ial,nsteps,nneig,ineig,
    c       Now grow from ial
            call growchain(n0,n,0,ial,iaf,lbackb,nneig,ineig,
    c       Save nn list for manipulations
              call trnsfi(ineiga(1,ia),ineig(1,ia),nneig(ia))
    c       Drop backbone atoms from nn list to avoid loopbacks
    c        do ia=n0,n
    c          write (77,6733) ia,(ineiga(i,ia),i=1,nneiga(ia))
    c6733      format(i5,' reduced nn=',30i4)
    c        end do
    c         If there is a sidechain, get its backbone too
    c             write (77,*) 'ia,in,ib,ibna=',ia,in,ib,ibna
                  call growchain(n0,n,iba,ibna,ial,nsteps,nneiga,ineiga,
    c             (First) monovalent neighbor
    c           Get the chain from ichainmax
    c           write (77,*) 'Longest sidechain:'
                call growchain(n0,n,iba,ichainmax,ial,nsteps,nneiga,ineiga,
    c         Print line(s)
    c             write (linep(ic0+1:ic0+18),1001) icis,
    c1001  format(' -',i4,' (',a4,a5,') ')
    c     write (77, *) 'GROWCHAIN n0,n,n00,iaf=',n0,n,n00,iaf
    c     write (77, *) 'GROWCHAIN nforbid,iaforbid=',nforbid,iaforbid
          call zeroiti(iparent,n0-1,n)
    c             Consider only atoms within the range
    c     nsteps is the number of vertices on the longest path, backtrack
    c      write (77,1000) n00,iaf,(ichain(ia),ia=1,nsteps+1)
    c1000  format(' n00,iaf=',2i4,(' chain:',25i4))
    c     Swap i-ia bond to the from n1 to n2 position
            call swapi4(ineig(n1,ia),ineig(n2,ia))
    c       ii=ineig(n1,ia)
    c       ineig(n1,ia)=ineig(n2,ia)
    c       ineig(n2,ia)=ii
            call swapi4(ibnd(n1,ia),ibnd(n2,ia))
    c       ii=ibnd(n1,ia)
    c       ibnd(n1,ia)=ibnd(n2,ia)
    c       ibnd(n2,ia)=ii
          character* 132 line(maxrec)
    c*****Sort atoms within a residue to follow an RTF order
          character*4 ires,reso,resn,atnam,sego,segn
          character*8 convdat
          character*200 sfilename
          common /savedat/ mxresdat,maxcondat,ifst(1000),ilst(1000),
          character*5 crdext
          common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
    c     print *,'sortatresrtf n=',n
          call trnsfi(indexo,index,n)
          call zeroiti(indexn,0,n)
    c         End of residue found
    c           Find the residue
                call findname(reso,ires,1,nres,ix,lench)
    c           print *,'nresproc,reso,ix=',nresproc,reso,ix
    c             Check the number of atoms
    c             Sort the residue
    c              write (77,7712) ifadb,iladb,
    c      -         (convdat(2,i)(1:4),i=ifadb,iladb)
    c7712          format('ifadb,iladb=',2i4,' convdat2=',/(10a5,/))
                    call leftadjust4(atnam,atnam)
    c               write (77,*) 'iaa,resnam,atnam=',iaa,reso,'*',atnam,'*'
    c6577  format(1x,a,' ifaslt,ilaslt=',2i5,' indexn:',(/,10i5))
    c                write (6,6577) 'Before compacting',ifaslt,ilaslt,
    c     -            (indexn(iaa),iaa=ifaslt,ilaslt)
    c               First compact the index list
                  call indexit(indexn,ifaslt,ilaslt,0)
    c     Now, indexn contains the new order (from 1 to n) or zero for
    c     database atoms that were found no match.
    c6711  format(1x,a,/,(30i4))
    c      write (6,6711) 'index',(index(i),i=1,n)
    c      write (6,6711) 'indexn',(indexn(i),i=1,n)
    c     Check for missed atoms
    c      write (6,6711) 'Cond indexn',(indexn(i),i=1,n)
    c       isegno should be the same within residues
    c      write (6,6711) 'Cond index',(index(i),i=1,n)
    c       Macromodel, rerrange the connectivity
    c           write (6,888) ia,ib,ic1,ic2,iold
    c888        format(' ia,ib=',2i3,' ic1,2=',2i3,' iold=',i3)
              call askstop(0)
          character*1 gcent,ans1
          character*4 ptype
          character*8 resn,resnam,atomn,atomnam
          character*4 ires
          character*8 convdat,grpinfo
          character*200 sfilename
          common /savedat/ mxresdat,maxcondat,ifst(1000),ilst(1000),
    c     print *,'PDBtommc resn=',resn,' atomn=',atomn
          char=0.0
          call leftadjustn(atomn,atomnam,8)
          call leftadjustn(resn,resnam,8)
    c     If residue name/atom name not found, check alternatives
    c       Original Amber
    c           Ask if oxy or deoxy NA
                call getname(ans1,len,
    c       Amber 94
    c       Charmm
    c???
    c     if (idigit(atomnam(1:1),1) .eq. 1) then
    c       atomnam(1:4)=atomnam(2:5)
    c       atomnam(5:5)=' '
    c     end if
    c       write (77,*) 'resnam,ires(i)=',resnam,ires(i),' if,l=',ifst(i),ilst(i)
    c            write (77,1077)
    c     -        resnam,atomnam,convdat(1,j)(1:4),convdat(2,j)(1:4)
    c1077        format(' rn,an=',a4,',',a4,'*',' cd1,2=',a4,',',a4,'*')
    c         If not found, try ACE, NME or DPOM or RPOM
    c         since Amber uses these groups as different residues
    c           Check if there are residue-name independent atomnames
          character* 132 line(maxrec)
          character*1 ans1
          character*80 inpline
          character*4 ires
          character*8 convdat
          character*200 sfilename
          common /savedat/ mxresdat,maxcondat,ifst(1000),ilst(1000),
          character*8 atnam
          character*200 filenames(5)
          character*57 RTFtyps
          common /explainRTF/ RTFtyps(4)
            call getname(sfilename,lsfilename,'Name of the conversion file',
    c       Read in conversion information
            call openfile(iuconv,0,'conversion information',22,'old',
    c         write (77,*) nconvdat,inpline
            call getrtfdat(nconvdat,5,iconvtyp,ntypedat)
    c      do i=1,nconvdat
    c        write (77,7711) i,(convdat(k,i),k=1,ncol)
    c7711    format(i5,5(2x,a8))
    c      end do
    c     Find the residue limits in convdat
    c         New residue found
    c     write (77,*) 'iresgen,nconvdat,nres=',iresgen,nconvdat,nres
    c     do ii=1,nres
    c       write (77,*) ii,ires(ii),ifst(ii),ilst(ii)
    c     end do
            call askyn('Do you want to de-regularize them',33,1,1,idereg,0,
              call fixrecform(line,index,n,3,inamcol1,inamcol2,
          character*80 linep,linech,lineprev
          character*4 ires,resnam,atnam,atnamo,potnam,attyp
          character*200 sfilename,rtffile,parmfile
          character*8 convdat
          common /savedat/ mxresdat,maxcondat,ifst(1000),ilst(1000),
          character*1 lc
          call openfile(15,0,'RTF',3,'old',rtffile,namlena,notfnd,0,1,1,0,0)
    c     First figure out if Amber, Charmm, or Gromacs
    c       Read Amber RTF
            call nextblank(linech,ic,80)
    c         write (77,*) 'resnam=',resnam
                call blankout(linech,1,80)
                call nextchar(linech,ic,80)
                  call nextblank(linech,ic,80)
                  call nextchar(linech,ic,80)
    c             write (77,*) 'atnam=',atnam
                      convdat(icl,nconvdat)='        '
    c               Residue name
                    convdat(1,nconvdat)(1:4)=resnam
    c               Atom name
                    convdat(2,nconvdat)(1:4)=atnam
    c               Atom type
                    call nextblank(linech,ic,80)
                    call nextchar(linech,ic,80)
                    convdat(3,nconvdat)(1:2)=linech(ic:ic+1)
                      call nextblank(linech,ic,80)
                      call nextchar(linech,ic,80)
                    call nextblank(linech,ic,80)
    c               Charge
                    convdat(4,nconvdat)(1:ic-ic1)=linech(ic1:ic-1)
    c               write (77,7711)nconvdat,(convdat(ic,nconvdat),ic=1,ncol)
    c7711            format(i5,7(1x,'|',a8,'|'))
    c         Search for DONE
    c       Read Charmm RTF
                  convdat(icol,nconvdat)='        '
    c           Residue name
                convdat(1,nconvdat)(1:4)=resnam
    c           Atom name
                call nextchar(linech,ic,80)
                convdat(2,nconvdat)(1:4)=linech(ic:ic+3)
    c           Atom type
                call nextchar(linech,ic,80)
                convdat(3,nconvdat)(1:4)=linech(ic:ic+3)
                  call readreal(linech,ic0,ic-1,charge)
                  charge=0.0
    c           print *,'RES ',resnam,' done'
                call nextchar(linech,ic,80)
                call nextblank(linech,ic,80)
    c       Read Gromacs RTF (.rtp) files
              call blankout(linech,1,80)
                  call nextblank(linech,ic,80)
                  call nextchar(linech,ic,80)
    c               Residue found
    c???            call getname4(ic,lineprev,resnam,80,1)
                    call getname4(ic,lineprev,resnam,80,1)
    c           Read record
    c             Residue ended
                  call getname4(ic,linech,atnam,80,1)
                  call getname4(ic,linech,potnam,80,-1)
    c               Residue ended
                      convdat(icol,nconvdat)='        '
    c               Residue name
                    convdat(1,nconvdat)(1:4)=resnam
    c               Atom name
                    convdat(2,nconvdat)(1:4)=atnam
    c               Atom type
                    convdat(3,nconvdat)(1:4)=potnam
                    call nextchar(linech,ic,80)
                    call nextblank(linech,ic,80)
                    call nextchar(linech,ic,80)
    c               print *,'names=',resnam,atnam,potnam
    c               print *,'charge=',charge,' ic=',ic,' lic=',linech(ic:ic)
          call askyn('Do you have an other RTF file',29,1,0,morefile,0,0)
            call openfile(16,0,'PARAM',5,'old',parmfile,namlenp,notfnd,
    c       Extract L-J eps and sigma values for atomtypes
    c         Amber
    c         Charmm
              call lastchar(linep,ic,80)
    c         Skip continuation line(s)
                call lastchar(linep,ic,80)
    c         Keep reading until blank line is found
                call nextchar(linep,icf,80)
    c             Read record
    c             Skip one number
                  call nextchar(linep,icf,80)
                  call nextblank(linep,icf,80)
    c             Read eps
                  call nextchar(linep,icf,80)
                  call nextblank(linep,icl,80)
                  call readreal(linep,icf,icl-1,epsilon)
    c             Read sigma
                  call nextchar(linep,icf,80)
                  call nextblank(linep,icl,80)
                  call readreal(linep,icf,icl-1,sigma)
    c             write (77,*) 'ntypedat,eps,sig=',
    c    -            ntypedat,eps(ntypedat),sig(ntypedat)
                call lastchar(linep,ic,80)
    c         Gromacs
            call askyn('Do you have an other PARAM file',31,1,0,morefile,0,
    c       write (77,2006) (attyp(i),eps(i),sig(i),i=1,ntypedat)
    c2006  format(' Atom type   Eps     sigma',/,(5x,a4,2x,f8.4,f8.4))
    c         Find type, extract eps & sigma
    c     write (77,*) 'nconvdat=',nconvdat
    c     do ncd=1,nconvdat
    c       write (77,7711) ncd,(convdat(ic,ncd),ic=1,ncol)
    c     end do
    c7711  format(i5,7(1x,'|',a8,'|'))
          character*(*)flag,form
          character*80 liner
          call blankout(liner,1,80)
            call lastchar(liner,lc,80)
          character*(*) name,list
    c     write (6,*) 'FINDNAME ifrst,lenlist,name=',ifrst,lenlist,name
    c      write (77,7711) (list(i),i=1,lenlist)
    c7711  format(10a5)
          character*(*) query,list(nlist)
          character*80 line
    c     print *,'PICKNAME nlist=',nlist
            call blankout(line,llist(i)+1,maxlen+1)
          call getint(query,lquery,0,1,nlist,ix,000)
          character*(*) namelist,label
    c     Checks list of names for duplicates
            call findname(namelist(i),namelist,1,i-1,ix,lnamelist)
          character*1 altcol(maxrec),inscol(maxrec)
          character*80 trtitle(32)
          character* 132 line(maxrec),blankline
          character*80 title
          character*4 segnames(mxres)
          character*8 resnamslv,atnames(maxrec),resnames(mxres)
          character*6 marker(16)
          character*5 crdext
          common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
          common /boxdat/ xtlabc(6),xtlabc0(6),box(3),box0(3),edge_gen(3,3),
          character*1 xyz
          common /axislab/ xyz(3)
          common /saveinfo/ nwrite,nfreeatwrite
          character*4 pflsv(100)
          character*8 namesv(100)
          character*4 charmmheader(2)
          character*11 charmmheadname(2)
    c     print *,'MAKETRAJ nslt,n,naslv=',nslt,n,naslv
    c     print *,'MAKETRAJ ioutrajtyp,ioppbc,icntrlr(11)=',
    c    -  ioutrajtyp,ioppbc,icntrlr(11)
    c       Create and write header
    c       Decide about cell info writing
    c         Charmm
                    call askyn('Do you want repeated box information',36,
    c               Input box size
                      call getreal(
                  call zeroitd(xtlabc,6)
    c         Amber
                call askyn('Do you want repeated box information',36,
    c           Use Charmm box size
    c           Input box size
                call getxyz('Initial box dimension in the ',29,
    c         Charmm
    c           call copybits(icntrlr(10),tsakma)
                call copyintgreal(icntrlr(10),tsakma,1)
                call getint('Number of fixed atoms to put in the header',
                call getint(
                call getint('Previous run steps to put in the header',39,
                call getint(
                call getint('Charmm version to put in the header',35,
                call getreal('Time step/fs to put in the header',33,tsfsdef,
    c           call copybits(tsakma,icntrl(10))
                call copyintgreal(icntrl(10),tsakma,0)
                call askyn(
    c           See if nfreeat has to be reduced when solvents are dropped
    c         Amber
    c         MMC
    c         if (naslv .eq. 3) then
    c           call getint('Number of solvent atoms per solvent molecule',
    c    -        44,999999,1,0,naslv,0)
    c         end if
              call extension(c,ih,0,1,nslt,cmin,cmax,shiftmmc,0,0,v)
    c       Create array with new order
                  c1(k,i-ndel)=c(k,indexs(i))
    c       Charmm
    c       print *,'MAKETRAJ icntrl(11),limic11=',icntrl(11),limic11
    c       Amber
    c       MMC
    c       Macromodel
            call writeconf(iout,inptp,iommod,inpcrdtyporg,nwrite,nwrite,
    c       Just write the shortened version (Xcluster)
          character*1 asterisk,altcol(maxrec),inscol(maxrec)
          character*4 segnames(maxrsd),segid4(maxrsd)
          character* 132 line(maxrec),blankline
          character*80 title,lineinp,trtitle(32)
          character*4 namin(maxrec),namout(maxrec),
          common /boxdat/ xtlabc(6),xtlabc0(6),box(3),box0(3),edge_gen(3,3),
          common /pbcrotmat/ torot_ac(3,3),torot_ca(3,3),tofac_ac,tofac_ca
          character*200 trajnam,outfile
          common /trajname/ trajnam,outfile,ltrajnam,namleno,ifirsttraj,
          common /columnlim/ incol(19),iidcol(19),iialtcol(19),iiinscol(19),
          character*5 crdext
          common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
          character*8 resnamslv,atnames(maxrec),resnames(maxrsd)
          character*6 marker(16)
          character*1 xyz
          common /axislab/ xyz(3)
          character*11 trajformatname
          common /trajectory/ nmmccheck,iftrajtyp(6),trajformatname(6)
          character*200 reffile
          character*1 ans,seqtyp,reftyp,chainid,chainid_prev
          character*4 resnam4,atnam4,resnam4n,atnam4n,segid,segid_prev
          character*8 resold,namold
    c     print *,'CONVTRAJ ioutrajtyp,inptrajtyp,inpcrdtyp,n=',
    c    -  ioutrajtyp,inptrajtyp,inpcrdtyp,n
          call indexit(indexs,1,n,0)
    c     Conversion - set up sorting
    c       Leave atom order unchanged
            call indexit(indexs,1,n,0)
              call select(line,nrecdel,idcol,asterisk,n,nslt,index,ixres,
    c       Read index file
            call quiz(seqtyp,ieqtyp,'s',' ',0,'sequence file type',18,0,5,6,
            call zeroiti(indexn,0,n)
    c         Open the .sort file
              call openfile(69,0,'index',5,'old',reffile,lreffile,
    c         Read new order from a PDB file' ATOM records
    c         Open the .pdb file
              call openfile(69,0,'index PDB',9,'old',reffile,lreffile,
    c         Open the Charmm .CRD file
              call openfile(69,0,'index CRD',9,'old',reffile,lreffile,
            close (69)
    c       Make a match
    c       Read output structure for atom and residue name info
            call openfile(69,0,'reference order',15,'old',reffile,lreffile,
              call quiz(reftyp,ireftyp,' ','template',8,
    c         Charmm .CRD
    c         PDB
              chainid_prev=' '
                  chainid=lineinp(iisegcol(1,ireftyp):iisegcol(2,ireftyp))
                    chainid_prev=chainid
            call getdupindex(nsegm,segid4,ixseg2)
            call initnamconv(noconv)
    c       Convert input names to output's convention
                call namconv(nrescol,resnam4,atnam4,resnam4n,atnam4n,nrch,
            call getdupindex(nsegm,segid4,ixseg1)
    c       cv and rprox will not be used in maketraj/writeconf calls
            call residue_contig(n,iresno,isegno,index1,ixseg1,indexo,
            call set_res_lim(iresno,n,ifres,ilres,nres,index1,1,maxrsd,
    c        do ir=1,nres
    c          do ia=ifres(ir),ilres(ir)
    c            write (77,8751) ir,ia,index1(ia),resin(index1(ia))
    c8751        format(' ir=',i4,' ia=',i6,' index1=',i6,' resin=',a)
    c          end do
    c        end do
            call residue_contig(nref,ires_ref,isegno_ref,index2,ixseg2,
            call set_res_lim(ires_ref,nref,ifres_ref,ilres_ref,nres_ref,
    c        do ir=1,nres_ref
    c          do ia=ifres_ref(ir),ilres_ref(ir)
    c            write (78,8752) ir,ia,index2(ia),resout(index2(ia))
    c8752        format(' ir=',i4,' ia=',i6,' index2=',i6,' resout=',a)
    c          end do
    c        end do
    c       Set up matching array
            call openfile(69,0,'order',5,'new',reffile,lreffile,notfnd,
            call zeroiti(indexn,0,n)
            call zeroiti(indexs,0,n)
    c         Find the next matching residue in the template
    c         print *,'Start checking ir,ir_ref=',ir,ir_ref
    c           if (ia .le. nref) write (6,1020) ia,resin(ia),namin(ia)
            close (69)
              call askyn('Do you want to try matching again',33,1,1,newm,0,
          call opentraj(c,1,inpt,inptrajtyp,n,ntitltr,trtitle,
              call askyn('Do you want to override the cell size/shape',43,
              call setpbccell('',0,edge,edge_gen,cell,ncell,cellalt,
              call trnsfr(cell0,cell,3*ncell)
    c       Generate cell info from trajectory cell sizes
    c           Skewed hexagon
              call trnsfr(edge,box0,3)
              call askyn(
    c           TO,  Amber to Charmm
              call askyn(
    c           TO, Charmm to Amber
              call trnsfr(edgexyz,edge,3)
    c       Recreate the cell from the first frame's cell size
            call crorgn(edgexyz,edge_gen,ioppbc,3,ncell,cell,cellalt,
            call trnsfr(cell0,cell,3*ncell)
            call askyn('Do you want to continue',23,1,-1,icont,20,0)
          call unitmat(trajrot)
            call askyn(
              call genrot(trajrot,pi,iax,angle)
              call askyn(
                call askyn('Do you want to select atoms for overlay',39,
                  call select(line,nrecdel,idcol,asterisk,n,nslt,index,
    c             indexdel: 0 for atoms to use for the overlay, 1 for the rest
    c             print *,'nfinalov=',nfinalov
    c             write (6,9173) (indexdel(i),i=1,nfinalov)
    c9173          format(80i1)
                  call masktolist(indexsup,indexdel,n,nfinalov,0)
    c             Get a condensed list of atoms selected to write
    c             print *,'nfinalov=',nfinalov,' nwrmax=',nwrmax,' nd=',ndel
              call trnsfr(co,c,3*n)
              call askyn(
              call setpbccell(
              call trnsfr(cell0,cell,3*ncell)
              call askyn('Do you want to try both TO orientations',39, -1,1,
            call matprod(trajrot,torot_ac,trajrot)
            call matprod(trajrot,torot_ca,trajrot)
            call setpbcdim(ioppbc,ixyzhex,ixyzexcld,ixyzincld,xyz)
            call setrepats(nmolslt,molsltlim,nslt,nneig,ineig,
            call quiz(ans,iansrun,'o',' ',0,
              call getxyz('',0,' component',10,999999.0,shiftpbc,0,0)
              call getint('Atom index of the atom to be at the center',42,
              call askyn(
            call prtcell(ioppbc,edge,edge_gen,r,vol,nw,1)
    c         Option to shift center to box corner
              call askyn(
    c         Option to shift center from box corner
              call askyn(
          call openfile(iout,0,'output trajectory',17,'new',outfile,namleno,
              call askyn('Do you want to use the structure title for both',
              call askyn('Do you want to add the structure title',
    c       Read a conformation
            call readtraj(inpt,inptrajtyp,mmctrajtyp,n,naslv,nwatr,nslt,
    c       r12=sqrt(dist2(c(1,molsltlim(3,1)),c(1,molsltlim(3,2))))
    c       if (r12 .gt. 40.0) then
    c         write (77,*) 'NMC=',nmc,' R12=',r12
    c         nr50=nr50+1
    c       end if
    c         See if this configuration is on the list
    c           List exhausted
                call askyn('Do you want to give an energy threshold (max)',
    c         Write a (converted) conformation
                call comparetop(c,n,nneig,ineig,iatnum,innlist,nslt,
    c           Reset into PBC
                call systemcenter(n,nmolslt,nmolsltnoion,molsltlim,c,c2,
    c           r12=sqrt(dist2(c(1,molsltlim(3,1)),c(1,molsltlim(3,2))))
    c           if (r12 .gt. 40.0) then
    c             write (78,*) 'NMC=',nmc,' R12=',r12
    c             nr50t=nr50t+1
    c           end if
    cDB            call debug_c(noutconf,line,index,nslt,cell,ncell,
    cDB     -        molsltlim(1,1),molsltlim(2,1),
    cDB     -        molsltlim(1,2),molsltlim(2,2),maxrec)
                call bestoverlay(nfinalov,indexsup,indexsup,co,c,atw,0.d0,
                call shiftmol(c,n,crm,c1,-1.0)
                call rotate_c(c1,n,rot,c2,'OVERLAY',7)
                call shiftmol(c2,n,c0,c,+1.0)
                call rotate_c(c,n,trajrot,c,'CONVTRAJ',8)
              call maketraj(noutconf,maxconf,c,c1,rprox,cv,ih,nslt,n,naslv,
    c     write (6,*) 'NR50=',nr50,' NR50t=',nr50t
    c1017  format(' Nconf=',i8,' rc1,rcmin1=',f7.1,f8.1,' rc2,rcmin2=',
    c     -  f7.1,f8.1,/,' n11-2,n21-2=',i4,3i5,' vtot1=',f10.1,' vtot2=',
    c     -  f10.2,' rmsd1=',f7.1,' rmsd2=',f7.1,' isw=',i1,/)
    c1020  format(' ERROR: no match for atom ',i5,' res=',a5,' name=',
    c     -  a5)
          common /nnwork/ cvsums(3,MAXREC),cvsums2(3,MAXREC),
    c     print *,'FILTERSLV NUMSOLV=',numsolv,' NSLT,N=',nslt,n,
    c    -  'ifilttyp=',ifilttyp
    c     print *,'FILTERSLV n,nslt,naslv=',n,nslt,naslv
          cang12min1=cos(ang12min1*3.141592/180.0)
          cang12min2=cos(ang12min2*3.141592/180.0)
          call zeroiti(ixdrop,0,numsolv)
    c     print *,'FILTERSLV numsolv,ifilttyp=',numsolv,ifilttyp
    c     print *,'FILTERSLV rsltmax,rcvmax=',rsltmax,rcvmax
            call calc_cv_rmin(c,ian,nslt,1,nslt,n,numsolv,naslv,iarepslv,
    c         write (78,8921) is,cv(is),rnearsq(is),ixdrop(is),rminmin
    c8921     format(i7,' cv=',f8.4,' rnsq=',f10.3,' ixdrop=',i2,
    c    -      ' rminmin=',f8.2)
              call calc_cv_rmin(c,ian,nslt,molsltlim(1,mol1),
                  call calc_cv_rmin(c,ian,nslt,molsltlim(1,mol2),
    c                   Solute anchor atoms are closer than r12max (15)
                        cang12min=cang12min1
    c                     Solute-protein line angle > ang12min
    c                       cv ratio is < cvrlim (3.5)
    c                       Calculate combined CV
                              cvsum(k)=cvsums(k,is)+cvsums2(k,is)
                            cv12=1.d0-
    c                         CV wrt both proteins is > cv12lim (0.6)
    c       Check for empty neighborhood
    c       do is=1,numsolv
    c         Calculate combined CV
    c         do k=1,3
    c           cvsum(k)=cvsums(k,is)+cvsums2(k,is)
    c         end do
    c         cv12s=1.d0-dsqrt(cvsum(1)**2+cvsum(2)**2+cvsum(3)**2)/
    c    -      dfloat(ncvsums(is)+ncvsums2(is))
    c         write (6,9681) is,ixdrop(is),cv12s,
    c    -      (c(k,nslt+(is-1)*naslv+iarepslv),k=1,3)
    c9681     format(i4,' ixdrop=',i2,' cv=',f8.5,' c=',3f10.5)
    c       end do
    c       Hydrogen-bond bridging water search
    c                 2nd Hbond was found to a different residue
    c     print *,'INDEX_UPDATE incr=',incr,' ix1=',index(incr+1)
    c     print *,'CALC_CV_ numsolv,naslv,iarepslv=',numsolv,naslv,iarepslv
    c     print *,'RSLTMAX,RCVMAX,CVMIN=',rsltmax,rcvmax,cvmin
    c     print *,'CALC_CV nslt1,nslt2,nslt=',nslt1,nslt2,nslt
          call cellpart(c,ian,itemp1,nslt1,nslt2,0.0,spacing,corner,ecell,
          call zeroit(cv,numsolv)
          call zeroiti(ianear,0,numsolv)
    c     do is=1,1
    c       Set the grid limits for search
            call trnsfr(cslv,c(1,incr+iarepslv),3)
    c       write (77,8791) is,ifar,cslv
    c8791   format(i7,' IFAR=',i2,' CSLV=',3f10.5)
    c       if (ifar .eq. 0) write (77,*) 'IXYZ1=',ixyz1
    c       if (ifar .eq. 0) write (77,*) 'IXYZ2=',ixyz2
    c       if (inside .eq. 3) then
    c         Check the box enclosinf the solvent
    c         ic=1+ix+iy*nxyz(1)+iz*nxyz(1)*nxyz(2)
    c         if (ifirst(ic) .gt. 0) then
    c           do ia=ifirst(ic),ilast(ic)
    c             ix=(c(1,indexs(ia))-xyzmin(1))/ecell(1)
    c             iy=(c(2,indexs(ia))-xyzmin(2))/ecell(2)
    c             iz=(c(3,indexs(ia))-xyzmin(3))/ecell(3)
    c             icc=1+ix+iy*nxyz(1)+iz*nxyz(1)*nxyz(2)
    c             write (77,*) 'IXYZ(c),ic=',ix,iy,iz,icc
    c           end do
    c         else
    c           write (77,*) 'EMPTY'
    c         end if
    c       end if
    c       Use all cells for debug
    c       do k=1,3
    c         ixyz1(k)=0
    c         ixyz2(k)=nxyz(k)-1
    c       end do
    c       ifar=0
            call zeroitd(cvsum,3)
    c                         write (77,*) 'd=',d
                              call norm(d,1.0)
    c                         write (77,*) 'dnorm=',d
                                cvsum(k)=cvsum(k)+d(k)
    c                         write (77,*) 'cvsum=',cvsum
            call trnsfrd(cvsums(1,is),cvsum,3)
    c       write (77,8945) is,sqrt(rnearsq(is)),cv(is),ncvsum
    c8945   format(i8,' rnear=',f8.2,' cv=',f6.3,' ncvsm=',i7)
    c     print *,'SET_RES_LIM n=',n
    c     print *,'RESIDUE_CONTIG n=',n
          call indexit(index,1,n,0)
          call mrgsrti(6,index,isr,n,ifa,ila,itmp1,itmp2,max)
          character*4 segid4(nsegm)
    c     print *,'GETDUPINDEX nsegm=',nsegm
          call indexit(index,1,nsegm,0)
          character*80 title,trtitle(32)
          common /boxdat/ xtlabc(6),xtlabc0(6),box(3),box0(3),edge_gen(3,3),
          character*200 trajnam,outfile
          common /trajname/ trajnam,outfile,ltrajnam,namleno,ifirsttraj,
          character*11 trajformatname
          common /trajectory/ nmmccheck,iftrajtyp(6),trajformatname(6)
          character*1 xyz
          common /axislab/ xyz(3)
    c     print *,'SPLITTRAJ inptrajtyp,mmctrajtyp=',
    c    -  inptrajtyp,mmctrajtyp
          call opentraj(c,1,inpt,inptrajtyp,n,ntitltr,trtitle,
          call askyn('Do you want to give the box size',32,0,1,nobox,0,0)
            call getxyz('Box size in the ',17,' direction (A)',14,
          call zeroiti(indexs,0,maxrec)
    c       Read a conformation
            call readtraj(inpt,inptrajtyp,mmctrajtyp,n,naslv,nwatr,nslt,
          call findlim(indexs,ifst,ilst,maxrec)
    c       Open output trajectory files
              call writeint(outfile,loutfile,it-1,lenw)
          call zeroiti(indexn,ifst-1,ilst)
    c       Read a conformation
            call readtraj(inpt,inptrajtyp,mmctrajtyp,n,naslv,nwatr,nslt,
            close (indexo(it))
          close (inpt)
          character*(*) label
          character*132 line(MAXREC)
          common /line_crd/ line,index
            call lastchar(line(index(ia)),lc,132)
            close(iu)
    cDB      common /DEBUG/ cx0(3,10000),
    cDB     -  cx1(3,10000),cx2(3,10000),cx3(3,10000),cx4(3,10000)
    c     print *,'SYSTEMCENT n,nmolslt,nslt,naslv=',n,nmolslt,nslt,naslv
    c     print *,'SYSTEMCENT nmolsltnoion,iacent,maxrsd,mxat=',
    c    -  nmolsltnoion,iacent,maxrsd,mxat
    c     print *,'SYSTEMCENT ncell,ixyzexcld,ixyzincld=',
    c    -  ncell,ixyzexcld,ixyzincld
    c     write (78,*) 'SYSTEMCENT'
    c     write (78,8822) (im,(molsltlim(k,im),k=1,3),im=1,nmolslt)
    c8822 format(' im=',i4,' molsltlim=',3i6)
    c     Bring together the solute molecules
    c     Find the center of the first solute molecule and shift it to <0,0,0>
    cDB      call trnsfr(cx0,c,3*nslt)
    c     print *,'IACENT,IMCENTER=',iacent,imcenter
    c           Default center: largest solute molecule's center
                call askyn(
                  call getint('Solute molecule number',22,1,1,imc_def,
              call extension(c,it1,0,molsltlim(1,imcenter),
    cxx       Calculate COM
              call zeroit(c00,3)
                  c00(k)=c00(k)+atw(ia)*c(k,ia)
                  c00(k)=c00(k)/atwsum
              call trnsfr(c00,c(1,molsltlim(3,imcenter)),3)
            call trnsfr(c00,c(1,iacent),3)
          call shiftmol(c,n,c00,c,-1.0)
    cDB      call trnsfr(cx1,c,3*nslt)
    c     call savepdb(88,'MOLEC_IMCENTER_CENTERED.pdb',27,c,n,0)
          call molreset(nmolfst,nmolslt,nmolshift,c,ct,molsltlim,it1,
    c     call savepdb(88,'AFTER_MOLRESET.pdb',18,c,n,0)
    cDB      call trnsfr(cx2,c,3*nslt)
    c     print *,'After  molreset'
    c     do is=nmolfst,nmolslt
    c       write (6,*) 'molsltlim=',molsltlim(3,is)
    c       call extension(c,it1,0,molsltlim(1,is),molsltlim(2,is),
    c    -      cmin,cmax,c0,1,0,v)
    c     end do
    c     Repeat, to bring in 2nd neighbor cell members
    c     call molreset(nmolfst,nmolslt,nmolshift,c,ct,molsltlim,it1,
    c    -  cell,ncell,cellalt,icellalt,imcenter,maxrsd,mxat)
    c     call extension(c,it1,0,1,nsltnoion,cmin,cmax,c0,0,0,v)
              call cofms(c,c0,nsltnoion,atw)
              call shiftmol(c,n,c0,c,-1.0)
              call arrsum(c0,c00,c0,3)
              call ionreset(nmolfst,nmolslt,nmolshift,c,molsltlim,it1,atw,
    cDB      call trnsfr(cx3,c,3*nslt)
    c     Check for ions outside the cell
            call pbcreset(c(1,molsltlim(1,is)),
              call pbcreset(c(1,molsltlim(1,is)),
    cDB      call trnsfr(cx4,c,3*nslt)
    c     Reset solvents
            call extension(c,it1,0,nslt+(iw-1)*naslv+1,nslt+iw*naslv,
    cx      call pbcreset(c(1,nslt+(iw-1)*naslv+1),naslv,c0,
    cx   -    cell,ncell,ixyzexcld,ixyzincld,img)
            call genimdist(c0,cell,1,ncell,img,d2)
    c     print *,'MOLRESET ncell=',ncell
          call extension(c,it1,0,molsltlim(1,imcenter),
    c     imgw=0
                call extension(c,it1,0,molsltlim(1,im),molsltlim(2,im),
                call combinevol(cmin,cmax,cmincent,cmaxcent,cmin12sav,
    c             write (6,8733) 1,cmin,cmax,cmincent,cmaxcent,cmin12sav,
    c    -          cmax12sav,volmin
    c           imgw=imgw+1
    c           call savepdb(89,'MOLRESET.pdb',12,c,molsltlim(2,nmolslt),
    c    -        imgw)
                  call trnsfr(ct,c,3*molsltlim(2,nmolslt))
                  call shiftmol(c(1,molsltlim(1,im)),
    c             imgw=imgw+1
    c             call savepdb(89,'MOLRESET.pdb',12,ct,molsltlim(2,nmolslt),
    c    -          imgw)
                  call extension(ct,it1,0,molsltlim(1,im),molsltlim(2,im),
                  call combinevol(cmin,cmax,cmincent,cmaxcent,cmin12,cmax12,
    c              write (6,8733) img,cmin,cmax,cmincent,cmaxcent,cmin12,
    c     -          cmax12,vol
    c8733          format(' IMG=',i2,' CMIN=',3f8.2,' CMAX=',3f8.2,/,
    c     -               ' CMINCENT=',3f8.2,' CMAXCENT=',3f8.2,/,
    c     -               ' CMIN12=',3f8.2,' CMAX12=',3f8.2,' VOL=',f12.2)
                    call trnsfr(cmin12sav,cmin12,3)
                    call trnsfr(cmax12sav,cmax12,3)
    c             write (6,9888) img,vol,imgopt,volmin
    c9888         format(' IMG=',i2,' V=',f12.3,' IMGOPT=',i2,
    c    -        ' VOLMN=',f12.3)
                  cmincent(k)=amin1(cmincent(k),cmin12sav(k))
                  cmaxcent(k)=amax1(cmaxcent(k),cmax12sav(k))
    c             Try the alternate (TO) cell orientation
                    call trnsfr(ct,c,3*molsltlim(2,nmolslt))
                    call shiftmol(c(1,molsltlim(1,im)),
                    call extension(ct,it1,0,molsltlim(1,im),molsltlim(2,im),
                    call combinevol(cmin,cmax,cmincent,cmaxcent,cmin12,
                      call trnsfr(cmin12sav,cmin12,3)
                      call trnsfr(cmax12sav,cmax12,3)
    c           print *,'IMGOPT=',imgopt
                    call shiftmol(c(1,molsltlim(1,im)),
                    call shiftmol(c(1,molsltlim(1,im)),
    c     close (89)
            cmin12(k)=amin1(cmin1(k),cmin2(k))
            cmax12(k)=amax1(cmax1(k),cmax2(k))
    c     print *,'MOLRESET ncell,ixyzexcld,ixyzincld=',
    c    -  ncell,ixyzexcld,ixyzincld
              call extension(c,it1,0,molsltlim(1,is),molsltlim(2,is),
              call cofms(c(1,molsltlim(1,is)),c0,
              call trnsfr(c0,c(1,molsltlim(3,is)),3)
            call genimdist123dim(c0,cell,1,ncell,ixyzexcld,ixyzincld,
    c         if (is .lt. 4) print *,'MOLRESET is=',is,' img=',img
              call shiftmol(c(1,molsltlim(1,is)),
    c         See if image is now in the central cell
              call shiftmol(c0,1,cell(1,img),c0,-1.0)
              call genimdist123dim(c0,cell,1,ncell,ixyzexcld,ixyzincld,
          character*1 ctyp
          character*41 question
    c     print *,'SETREPATS  maxneig,maxrsd=', maxneig,maxrsd
    c     Establish solute molecule centers (if required)
          call zeroiti(moltyp,0,3)
            call askyn(
    c     molsltlim(1,im), molslt(2,im): first and last atom of solute molecule im
    c     molsltlim(3,im): atom number of topological center;
    c                    zero if geom center was selected; -1 if COM is selected
    c     Optionally establish representative atoms
          call quiz(ctyp,ictyp,'g',' ',0,'molecular center',16,0,5,6,25)
              call findtcent(ineig,nneig,it1,it2,molsltlim(1,is),
    c#    MMC routine 463 lstmod: 01/20/05
    c*****Find the topology center of a molecule
    c     print *,'FINDTCENT n0,n,maxneig,maxat=',n0,n,maxneig,maxat
    c         First grow neighbour list from n0
    c         Now grow from list(ilist)
    c               Consider only atoms within the range
    c       Now list(ilist) is one of the farthest from ic
    c       write (77,1000) n0,n,ir,list(ilist),nsteps
    c1000  format(' FINDCENT: Center search for atoms ',i5,' - ',i5,
    c     -  ': run',i2,' End point: ',i5,' Steps: ',i4)
    c     nsteps is the number of vertices on the longest path, backtrack
    c     nsteps/2 to get the center
          character*4 pflsv(100)
          character*8 namesv,resnamslv
          character*5 crdext
          character*80 question
          common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
    c     print *,'UPDATESOLVENTS nmax,n,naslv=',nmax,n,naslv
            call askyn('Is the solvent water (O,H,H)',28,1,1,iwat,0,0)
              call getreal('Oxygen charge',13,-0.834,qox,0,0)
                call getname(pflsv(1),len,question,31,4,'',0,0,0,0)
                call getname(pflsv(2),len,question,32,4,'',0,0,0,0)
              call getint('Number of atoms in a solvent molecule',37,
                call getint(question,33,999999,1,99,iasv(i),00)
                  call getname(namesv(i),len,question,29,5,'',0,0,0,0)
                  call getreal('Charge',6,999999.0,qsv(i),0,0)
                  call getname(pflsv(i),len,question,35,4,'',0,0,0,0)
              call askstop(0)
          character* 132 line(maxrec),blankline
          character*1 altcol(maxrec),inscol(maxrec)
          character*4 pflsv(100),segnames(maxrsd)
          character*8 resnamslv,atnames(maxrec),resnames(maxrsd),namesv(100)
          character*6 marker(16)
          character*80 trtitle(32),title
          character*200 inpfile
          common /logging/ logfile,ipredict
          character*5 crdext
          common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
          character*200 trajnam,outfile
          common /trajname/ trajnam,outfile,ltrajnam,namleno,ifirsttraj,
          character*11 trajformatname
          common /trajectory/ nmmccheck,iftrajtyp(6),trajformatname(6)
    c     print *,'UNPACKTRAJ inptrajtyp,mmctrajtyp=',inptrajtyp,mmctrajtyp
          call opentraj(c,1,inpt,inptrajtyp,n,ntitltr,trtitle,
            call getname(inpfile,namleni,
            call getname(inpfile,namleni,'Output coordinate file',22,200,
    c     write (6,7866) numsel,(iconfsel(i),i=1,numsel)
    c7866 format(' NUMSEL=',i4,' ICONFSEL=',20i6)
    c       Read a conformation
            call readtraj(inpt,inptrajtyp,mmctrajtyp,n,naslv,nwatr,nslt,
                call comparetop(c,n,nneig,ineig,iatnum,innlist,n,1,
    c           List-directed unpacking - check list
    c             Selection found - increment nextconfsel
    c           Full unpacking
                  call askyn(
    c           Write a conformation
                    call filenamenum(inpfile,namleni,outfile,nl1,ifnumw,+2)
                  call openfile(20,0,'output',6,'new',outfile,nl1,
                    call askyn(
                    call openfile(20,0,'output',6,'new',outfile,nl1,
    c             Add extra solvents
                  call updatesolvents(iaskatnum,1,nmax,n,naslv,iasv,namesv,
                call writeconf(20,inpcrdtyp,iotyp,inpcrdtyporg,nmax0,n,nslt,
          close(inpt)
          character*(*) prompt
          character*80 trtitle(32)
          character*200 inpfile
          character*4 chhd
          character*11 trajformatname
          common /trajectory/ nmmccheck,iftrajtyp(6),trajformatname(6)
          common /boxdat/ xtlabc(6),xtlabc0(6),box(3),box0(3),edge_gen(3,3),
          character*5 crdext
          common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
          character*12 recordn
    c     print *,'OPENTRAJ iverb,notitprint,irept,iopen_rep=',
    c    -  iverb,notitprint,irept,iopen_rep
    c     print *,'OPENTRAJ inptrajtyp,mmctrajtyp=',inptrajtyp,mmctrajtyp,
    c    -  ' ifirst,ilast=',ifirst,ilast
    c     print *,'OPENTRAJ icont_traj=',icont_traj,' IOUT_LIM=',iout_lim
          call blankout(recordn,1,12)
            call askyn('Do you want to change it',24,1,-1,newn,0,0)
    c       Open a Charmm trajectory
    c       print *,'Number of title lines=',ntitltr
                call blankout(trtitle(i),1,ntitltr)
    c             Check for NAMD origin
              call askstop(1)
              call askstop(0)
    c       Free atom array
    c         Test for presence of cell information
    c       Rewind and skip back (to be on the safe side)
    c         Cell information was found
              call trnsfrd(xtlabc0,xtlabc,6)
    c           Skewed hexagon
    c           Test for cell data before 2nd config
    c       Reposition trajectory
    c       Open an Amber trajectory
              call askyn('Do you want to change it (to include solvents)',
            call blankout(trtitle(1),1,80)
    c       Find first box info record
            call zeroit(z,10)
            call countzeros(z,10,nzero1)
    c       Check if box info was present after the 2nd config.
    c         Save the first box and skip to the end of the second frame
              call trnsfr(box0,z,3)
    c         See if there is second box info
              call zeroit(zz,10)
              call countzeros(zz,10,nzero2)
    c       Reposition trajectory
    c         MMC trajectory
                call binhst_type(inpt,c,ibox,istuner,ieof,6,maxrec)
                  call askstop(1)
                  call askstop(1)
            call blankout(trtitle(1),1,80)
    c       Amber CDF
    c           Set PBC to cubic, ask for size
                call pbcsize(ioppbc,edg,npbc)
    c         Recreate the cell from the first frame's cell size
              call crorgn(edg,edge_gen,ioppbc,3,ncell,cell,cellalt,
              call trnsfr(cell0,cell,3*ncell)
              call trnsfr(cell0,cell,3*ncell)
    c           See if a list is given for snapshots to be read
                call askyn('Do you have a list of configurations to read',
                  call getlist(iconfsel,numsel,1,999999,1,maxconfsel)
    c     Read header and coordinate record from a binary history file
    c     print *,'BINHST_READ maxat=',maxat
    c2000  format(' nwat,nats=',2i11,' nmc=',i10,' ia0,ia1=',2i6,
    c     -  ' nsv=',i2,/,
    c     -  ' etot=',e13.5,' tesi=',4e12.5,/,' cplpar=',e12.5)
    c     Determine if MMC binary trajectory file has box and/or tuning info
          character*4 chhd
    c     print *,'BINHST_TYPE maxat=',maxat
    c     Check for mistaken Charrm dcd
          call binhst_read(ihist,c,etoto,nmc,nwatr,ieof,0,maxat)
    c     New header was found
    c     Tuning and/or box record found
          call binhst_read(ihist,c,etoto,nmc,nwatr,ieof,0,maxat)
    c     Tuning record was found - no box info
    c     There must be box information
          call binhst_read(ihist,c,etoto,nmc,nwatr,ieof,0,maxat)
    c     No tuning info
    c     There is be tuning iformation too information
          character*11 trajformatname
          character*80 trtitle(32)
          character* 132 line
          character* 200 trajnam,trajnam_n(2)
          common /boxdat/ xtlabc(6),xtlabc0(6),box(3),box0(3),edge_gen(3,3),
          character*1 separatorchar
          common /filenuminfo/ iaskunderscore,separatorchar
          character*1 xyz
          common /axislab/ xyz(3)
    c     print *,'READTRAJ inpt=',inpt,' nslt,n=',nslt,n,
    c    -  ' inptrajtyp=',inptrajtyp
    c       Just use input edge info
    c         Binary MMC
    c           Skip tuning info
    c           ASCII MMC
    c           Annotated ASCII MMC
    c         PDB
                call blankout(line,1,80)
                  call checkforetot(6,line,ninconf,etot,ietotread,iverbconf)
    c         Charmm CRD
              call blankout(line,1,80)
                call checkforetot(1,line,ninconf,etot,ietotread,iverbconf)
    c       Amber CDF
    c       Just in case, generate Charmm edges
    c       print *,'n,nslt=',n,nslt
    c       Generate the next trajectory name
            call nextnames(trajnam,ltrajnam,trajnam_n,ltrajnam_n,ntraj)
              call opentraj(c,0,inpt,inptrajtyp,n,ntitltr,trtitle,inpcrdtyp,
          character*200 trajnam,trajnam_n(2)
          character*80 ext
    c     Find out run number/version number category and generate new name options
    c     Number read
    c       x.nr.ext
            call writeint(trajnam_n(2),lr,nr,len)
    c       x.$_nv.ext
    c       x.nr_nv.ext
            call writeint(trajnam_n(1),lc,nv,len)
            call writeint(trajnam_n(2),lr,nr,len)
            call writeint(trajnam_n(2),lr,nv0,len)
    c     x_nv.ext
          call writeint(trajnam_n(1),lc,nv,len)
          call writeint(trajnam_n(2),lc,nv,len)
    c     x.ext
          character* 132 line
    c     print *,'DB ECHECK lmarker=',lmarker,' line=',line(1:40)
          call nextchar(line,icol,80)
          call nextblank(line,icol,80)
            call nextchar(line,icol,80)
            call nextblank(line,icol,80)
          common /connatdat/ ramax(99),ramax2(99),hlimfac,ianfg(99),
    c     Check if list of neighbors actually represent bonded pairs
    c     print *,'COMPARETOP n,nslt,nslv,ncell=',n,nslt,nslv,ncell
    c       See how many solvents are outside the cell
              call genimdist(c(1,nslt+(iw-1)*nslv+1),cell,1,ncell,img,d2)
            call askstop(1)
          character*200 inpfile,outfiletmp
          character* 132 line(maxrec),ansline
          character*1 aanames1
          character*2 mmodtoamb
          character*3 aanames3
          common /atnamcon/ mmodtoamb(100),aanames1(58),aanames3(58),
          character*2 iatnm2
          common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
          character*5 crdext
          common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
          character*1 HMtyp,gcent,hml(3),s1
          character*2 ambtyp
          character*4 potnam
          character*8 resnam,atomnam
          character*80 prtline
    c     Generate output filenames
          call changeext(inpfile,outfiletmp,namleni,namlentmp,'onm',3,0,0)
          call openfile(21,0,' ',1,'new',outfiletmp,namlentmp,notfnd,
          call quiz(HMtyp,iHMtyp,' ','ONIOM',5,'region definition mode',22,
            call getlist(indexn,nhigh,1,nslt,1,nslt)
            call getlist(indexo,nmid,1,nslt,1,nslt)
            call getint('Atom number at the center of the quantum part',45,
            call getreal('Radius of the H region',22,5.0,rh,1,56)
            call getreal('Radius of the M region',22,10.0,rm,1,56)
            call askyn(ansline(1:37),37,1,-1,ihml(i),0,0)
            call zeroiti(mmtype,0,nats)
            call askyn('Do you want to prevent non C-C bond break',41,1,-1,
              call nnlist(nslt,nslt,naslv,n,iatnum,ifchrg,c,nneig,nneiga,
    c                 Add ja to the inner list
    c                 Add ja to the inner list
          call trnsfi(indexs,indexn,nhigh)
          call trnsfi(indexs(nhigh+1),indexo(nhigh+1),nmid)
    c     indexn,indexo,indexs contain the list of H, and M, atoms, resp.
    c       indexs contains the atomindices sorted by H, M, L
            call indexit(indexs,1,nats,0)
          call askyn('Do you want to reorder by ONIOM type',36,1,-1,iord,0,
          call zeroit(qhmlsum,3)
                call readint(line(index(ia)),ipotcol1,ipotcol2,mmodtyp,4,1,
              call PDBtommc(resnam,atomnam,potnam,qslt,ia,gcent,igr,
    c_C-CT-q   3f13.7,2x,a1
    c          ^ col 21
          call askyn('Do you want to generate topology input too',42,
            call nnlist(nslt,islvw,naslv,n,iatnum,ifchrg,c,nneig,nneiga,
            call bondord(iatnum,mmtype,n,nneig,ineig,nhneig,ibnd,maxng,c,
              call askyn('Do you want to add H-H bonds to the solvent',43,
          close(21)
          call trnsfi(isegno,ih,n)
          character* 132 line(maxrec)
          character*2 iatnm2
          common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
          character*5 crdext
          common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
    c     Obtain charges (if available)
              call readreal(line(index(ia)),iqcol1,iqcol2,qslt)
            call getint('The charge on the molecule',26,999999,1,0,iqsum,0)
    c     Obtain SM1A atomtypes
    c     Write header
               call writeline(iout,line(i),7,86,0)
               call writeline(iout,line(i),87,132,0)
    c     Write coordinates
          close (iout)
          character*1 asterisk,dssplab(maxrsd),altcol(maxrec),inscol(maxrec)
          character*4 segid4(nsegslt),segnames(maxrsd),extnam1,extnam2,
          character*8 resnam,rn,resnamslv,atnames(maxrec),resnames(maxrsd)
          character*6 marker(16)
          character*8 version
          character*24 askcolcode
          character*80 label2d(mx2d)
          character*200 inpfile,analfile,analfile1,analfile2,analfile3,
          character* 132 line(maxrec),blankline
          character*2 iatnm2
          common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
          character*11 trajformatname
          common /trajectory/ nmmccheck,iftrajtyp(6),trajformatname(6)
          character*200 trajnam,trajnam2
          common /trajname/ trajnam,trajnam2,ltrajnam,ltrajnam2,ifirsttraj,
          character*80 title,trtitle(32)
          character*4 namfcg
          character*4 tanames
          character*8 tnames
          common /tordat/ ntorn,tanames(4,28),tnames(28)
          common /connatdat/ ramax(99),ramax2(99),hlimfac,ianfg(99),
          common /boxdat/ xtlabc(6),xtlabc0(6),box(3),box0(3),edge_gen(3,3),
          common /colorinfo/ ncolcode,maxcolcode
          common /graphics/ npixx,npixy,maxpixx,maxpixy,idwmain,idwplot,
          common /rotmat/ matrot0(4,4),matrot(4,4),nomat0
          character*1 xyz
          common /axislab/ xyz(3)
          common /columnlim/ incol(19),iidcol(19),iialtcol(19),iiinscol(19),
          character*5 crdext
          common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
          character*4 pflsv(100)
          character*8 namesv(100)
          character*9 decidebend
          character*27 qq
    c     Maximum number of bridge anchor atoms: MAXBRDIGEATOM
    c     Maximum number of bridge destination atoms: MAXBRDIGETYPE
    c     Maximum number of H bonds in a bridge:  MAXBRIDGELEN
          common /bridges/ ianchor(MAXBRIDGEATOM),
    c     Maximum number of atom pairs to calculate the distance distribution
          character*1 typc
          character*21 ssname
          common /dsspnames/ lssname(9),ssname(9),typc(9)
    c     Pseudorotation calculation
    c     All arrays in prokink are of length MAXHX
          common /prokink/ icab(MAXHX),icaa(MAXHX),icb(MAXHX),ica(MAXHX),
          common /analparm/ nsltref_f,nsltref_l,rcut_cv,icvtyp
    c     All arrays are of length maxframe=MAXFRAMES > MAX2D !!
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
          common /logging/ logfile,ipredict
          character*1 separatorchar
          common /filenuminfo/ iaskunderscore,separatorchar
          character*2 ap_pa,in_ex
          common /signs/ itmem,normhx,isg2(MAXNHX2),memdir(MAXNHX),ap_pa(3),
    c     common /hxtrack/ ihxt1,ihxt2,torsav1(3,4),torsav2(3,4),x1(3),x2(3)
    c     Array of lenght maxconfsel (set to MAXFRAMES)
    c     Helix axis analysis arrays
    c     Arrays of length involving MAXNHX
          character*24 hxhxlab(MAXNHX2)
          character*1 axdirchar(MAXNHX)
          character*1 analtyp,analinp,psrtyp,ansrun,marks(9),unit,ans
          character*1 tmchar,tmchars(2)
          character*4 ext1(45),ext2(45),ext3(45),atnam,bondlab(5)
          character*8 brslv
          character*6 hxoklab(3),crdexti
          character*7 mark0,mark1,mark2
          character*11 xtrajlabs(4),xtrajlab
          character*20 helixcklab
          character*22 bondname(6)
          character*30 shiftlab
          character*80 question,linein,pstitle,system
          character*34 distlab,distlab_cc
          character*36 restitle
          character*22 atomdist
          character*25 prokinklab(5),helixang(6),helixrlab(14),rmsdlab(4),
          character*29 resrange
          character*30 talab(MAXCOPY1)
          character*80 plotdescr
          character*100 hostname
          character*200 trajnam1,trajnamr1,trajnamr2
    c     ianaltyp=1:  (S) Neighbor, bond, angle and torsion list
    c     ianaltyp=2:  (S) 1-4 statistics
    c     ianaltyp=3:  (S) Functional group and backbone list
    c     ianaltyp=4:  (S) Bond length statistics
    c     ianaltyp=5:  (S) Hydrogen-bond list
    c     ianaltyp=6:  (S) Hydrophobic bond list
    c     ianaltyp=7:  (S) Salt bridge list
    c     ianaltyp=8:  (S) Calculate residue distances
    c     ianaltyp=9:  (S) Calculate a PBC-adjusted distance
    c     ianaltyp=10: (S) Check for potentially unphysical contacts
    c     ianaltyp=11: (S) Pseudorotation angle calculation
    c     ianaltyp=12: (S) Calculate Proline kinks
    c     ianaltyp=13: (S) Hydropathy labeling
    c     ianaltyp=14: (S) Circular variance labeling of solute and solvent
    c     ianaltyp=15: (S) Circular variance residue-residue plot
    c     ianaltyp=16: (S) DSSP secondary structure assignment
    c     ianaltyp=17: (S) Hydrogen-bond bridge analysis
    c     ianaltyp=18: (S) Ramachandran plot
    c     ianaltyp=19: (S) Torsion dial plots
    c     ianaltyp=20: (S) Delphi map annotation
    c     ianaltyp=21: (S) Helix axis directions
    c     ianaltyp=22: (T) 1-D RMSD and residue RMS fluctuations
    c     ianaltyp=23: (T) 2-D RMSD
    c     ianaltyp=24: (T) Cross RMSD
    c     ianaltyp=25: (T) Residue correlation matrix calculation
    c     ianaltyp=26: (S) Atom-atom distance (distribution) calculation
    c     ianaltyp=27: (S) Solvation shell volume calculation
    c     ianaltyp=28: (S) Principal axis calculation
    c     ianaltyp=29: (S) Radius and dipole calculation
    c     ianaltyp=30: (S) Summarize Amber energy partition table
    c     ianaltyp=31: (S) Adjacency matrix analysis
    c     ianaltyp=32: (S) Angle dial plots
    c     ianaltyp=33: (S) Molecule-molecule distance list
    c     ianaltyp=34: (T) Heavy atom contact list
    c     ianaltyp=35: (S) Calculate eigenvectors from input matrix
    c     ianaltyp=36: (S) Compare residue-residue average distance matrices
    c     ianaltyp=37: (S) Compare residue-residue bond matrices
    c     ianaltyp=38: (S) Compare residue RMSF values
    c     ianaltyp=39: (T) Atom-atom SD from trajectory
    c     ianaltyp=40: (T) Solvent filtering
    c     ianaltyp=41: (T) Mutually proximal contact list
          call testconst(0,1,2,0.0,1.0,2.0,6,nfail,1,'ANST')
    c     do ir=1,nresslt
    c       write (77,*) 'ir=',ir,' iresf,l=',ifres(ir),ilres(ir)
    c       do ia=ifres(ir),ilres(ir)
    c         write (77,*) line(index(ia))(1:78)
    c       end do
    c     end do
    c     stop
          call indexit(ixshuffle,1,MAX2D,0)
          call indexit(ixshuffleref,1,MAX2D,0)
    c     print *,'numsel,maxconfsel=',numsel,maxconfsel
    c     print *,'maxrepconf,maxng,maxbox=',maxrepconf,maxng,maxbox
    c     print *,'maxrec,maxrsd=',maxrec,maxrsd
          call trnsfr(cres,c,3*n)
          call nextchar(title,ifc,1000)
          call lastchar(title,ltitle,80)
          call blankout(question,1,80)
            call getatnumlist(n,iatnum,ifchrg,ialist,icatlist,ixlist,nanos)
    c       Get residue name list
              call quiz(analtyp,ityp,' ',' ',0,
              call quiz(analtyp,ityp,' ',' ',0,'bond tracking',13,
              call quiz(analtyp,ityp,' ',' ',0,
              call quiz(analtyp,ityp,' ',' ',0,
              call quiz(analtyp,ityp,' ',' ',0,'RMSD calculation',16,
              call quiz(analtyp,ityp,' ',' ',0,'distance analysis',17,
                call askyn('Do you have an input covariance matrix',38,1,-1,
            call blankout(resrange,1,29)
              call quiz(ansrun,nntyp,' ',' ',0,'Bond information source',23,
                call top_to_bond(nntyp,nneig,nhneig,ineig,iatnum,n,0,
                call askyn('Do you want to change bond thresholds',37,1,
    c           Update bond thresholds
                    call getreal(linein,23,ramax(ian),ramaxnew,1,000)
    c         These options do not work on trajectories
    c         These options only work on trajectories
              call askyn('Do you want to analyze a trajectory',35,1,
                call quiz(unit,iframeunit,'f',' ',0,'trajectory unit',15,
                  call getreal(linein,45,1.0,framefac,1,113)
    c         icharges=1: charge in iqcol1:iqcol2; icharges=2: charges read
                call readcharges(nread,nslt,n,charge,iatnum,isv,
                  call checkreschargesum(nslt,iresno,isegno,line,index,
                  call askstop(1)
            call zeroiti(nhbdist,0,mxbonds)
            call zeroit(rhbdist,mxbonds)
            call indexit(indexs,1,MAX2D,0)
            call indexit(index2d,1,MAX2D,0)
            call indexit(ixclst,1,MAX2D,0)
                call askyn('Do you want to write a grid potential file',42,
                call askstop(0)
    c           Extend arrays with more solvents than the input has
    c           print *,'nadd=',nadd,' iasv=',iatnum(nslt+1),iatnum(nslt+2),
    c    -         iatnum(nslt+3)
                    charge(nslt+1)=-0.834
                    charge(nslt+2)=0.417
                    charge(nslt+3)=0.417
                    call getint(linein,45,999999,1,0,naslv,0)
                      call getint(linein,32,999999,1,99,iatnum(nslt+ia),0)
                    charge(n+(is-1)*naslv+ia)=charge(nslt+ia)
                call askyn(
    c         Open output file(s)
              call strip_cext(analfile,namleni,namleno,lenext)
    c         If output is a pdb file, add .pdb extension
    c           Select circular variance type
                call quiz(ansrun,icvtyp,'o',' ',0,
              call openfile(iw0,0,'analysis',8,'new',analfile,namleno,
              call datprt(iw0,version,1,mark0,lmark0,hostname,lhostname,
    c             Second file is PS - append the .ps to the name
    c             Second file is PDB - append the .pdb to the name
    c             Replace extension
            call datprt(6,version,1,' ',1,hostname,lhostname,iheadnode,1)
              call askyn(
                call askyn(
              call askyn('Do you want to draw arcs in the dial plots',42,1,
              call askyn('Do you want to print bond angles',32,1,-1,iangpr,
              call askyn('Do you want to print torsion angles',35,1,-1,
              call findchiral(nslt,iatnum,nneig,nhneig,ineig,indexa,maxng)
                call openps(iw1,xm,ym,title,ltitle76,'Ramachandran plot',
                  call askyn('Do you want dial plots for all residues',39,
                    call indexit(ixselres,1,nres,0)
                  call getint('Number of residues to track',27,0,1,maxrp,
                        call getint(question,38,isegdef,1,nsegm,isegix,00)
                      call getint(question,38,999999,1,nres,iresix,00)
                      call findsegres(isegno,iresno,ixres,1,nslt,isegix,
                        call leftadjust4(atnam,atnam)
                      call ca_to_bb(icafound,iresno,nneig,ineig,index,line,
                  call openfile(iw2,0,'dial plot',9,'new',analfile2,
                  call datprt(iw0,version,1,mark0,lmark0,hostname,lhostname,
                  call askyn(
                call ramachandran_init(n,ixres)
    c           Torsion angles
                call torslistinp(ixtor1234,talab,ltalab,ntorsel,inpcrdtyp,
    c           Bond angles
                call getint('Number of angles to track',25,0,1,49,
                  call getintline(question,35,1,nslt,ixtor1234(1,it),3,0)
    c             print *,' ixtor:', (ixtor1234(k,it),k=1,4)
                    call askyn('Do you want to use this angle',29,1,-1,iok,
                  call blankout(talab(it),1,30)
                call getint('Number of dials to draw in a line',33,
                call quiz(ans,itorcorr,'j',' ',0,
    c         Hydrogen-bond bridge analysis
              call set_hbondlims(hbf0,hblimfac,angm0,angmin,iw0)
    c             Check if residue exists, find # of atoms
              call getint('Maximum number of bridge residues in a bridge',
              call gethbanchordef(line,index,nslt,ixres,iresno,iatnum,
              call zeroiti(nbridgetype,0,nanchor*MAXBRIDGELEN)
              call askyn('Do you want to print the bridges',32,1,ipb,
              call getint('Minimum length of bridge to list',32,1,1,
                call getint('Minimum percent of frames present to list',41,
              call zeroiti(lpath,0,MAXBRIDGETYPE*MAXBRIDGELEN*
                call askstop(1)
                call askyn(linein,llinein,1,-1,ibondcorr,0,0)
                call askyn(linein,llinein,1,-1,iresbondcorr,0,0)
                call askyn(linein,llinein,1,-1,ibondprint,0,0)
                call askyn(
                    call getname(trackfile,ltrackfile,'Name of track file',
                    call openfile(iout_track,0,'previously written track',
                    call readtrack(iout_track,iw0,30,nbfound,nbresfound,
                call getint(
    c           Set up H bond analysis
                call set_hbondlims(hbf0,hblimfac,angm0,angmin,iw0)
                  call askyn('Do you want to include the solvents',35,1,-1,
                  call gethbanchordef(line,index,nslt,ixres,iresno,iatnum,
                  call indexit(it3,1,MAXBONDS,0)
                  call indexit(it4,1,MAXBONDS,0)
    c           Set up hydrophobic contact analysis
                call getreal('Hydrophopbic bond length limit',30,5.0,
                call getint('Minimum number of carbon-bonded hydrogens',41,
                call gethphanchordef(line,index,nslt,iresno,iatnum,charge,
                call extend_nnlist(nneig,ineig,npneig,nslt,maxng,maxrec)
    c            do ia=1,nslt
    c              atnam=line(index(ia))(inamcol1:inamcol1+3)
    c              write (77,7719) ia,atnam,nneig(ia),npneig(ia),
    c     -          (ineig(in,ia),in=1,npneig(ia))
    c            end do
    c7719        format(i5,1x,a,' nn=',i3,' nnp=',i6,' in=',(10i5))
    c           Set up salt-bridge contact analysis
                call getreal('Salt-bridge length limit',24,5.0,rsltbmax,1,2)
                call getsltbanchordef(line,index,nslt,iresno,iatnum,charge,
                call extend_nnlist(nneig,ineig,npneig,nslt,maxng,maxrec)
    c           Set up heavy-atom VdW contact analysis
                call getreal('Heavy-atom distance threshold',29,5.0,
                call gethphanchordef(line,index,nslt,iresno,iatnum,charge,
                call extend_nnlist(nneig,ineig,npneig,nslt,maxng,maxrec)
    c            do ia=1,nslt
    c              atnam=line(index(ia))(inamcol1:inamcol1+3)
    c              write (77,7719) ia,atnam,nneig(ia),npneig(ia),
    c     -          (ineig(in,ia),in=1,npneig(ia))
    c            end do
    c7719        format(i5,1x,a,' nn=',i3,' nnp=',i6,' in=',(10i5))
                  call askyn('Do you want to print the heavy-atom contacts',
    c           Set up mutually proximal contact analysis
                call getmpxbdef(nslt,indexa,indexov,indexn,segid4,iresno,
                call getreal('Maximum distance for contact',28,99999.0,
    c           call extend_nnlist(nneig,ineig,npneig,nslt,maxng,maxrec)
    c         Set up residue contact list calculation
              call modrepats
              call getreal('Threshold distance with representative atoms',
              call getreal('Threshold distance with closest approach',40,
              call askyn(
              call getresrange(nsegm,indexs,isegno,ixres,iresno,ifres,
              call getresrange(nsegm,indexs,isegno,ixres,iresno,ifres,
                  call askyn(
                  call openfile(iw4,0,'average distance matrix',23,'new',
                  call openps(iw4,500.0,500.0,' ',1,' ',1,inpfile,0,
    c           The ranges are disjoint
                call openfile(iw2,0,'contact list',12,'new',analfile2,
                call datprt(iw2,version,1,mark0,lmark0,hostname,lhostname,
                call askyn(
                  call askyn('Do you want to include solvent neighbors',40,
              call zeroiti(irescount1,0,nresslt)
              call zeroiti(irescount2,0,nresslt)
              call zeroiti(irescount3,0,nresslt)
    c         Set up molecule-molecule distance calculation
              call askyn(
    c           Compare two residue distance matrices
                call compare_rrdist(resnames,nrescol,itemp1,temp,irefres1,
                call askyn(
    c           Compare two residue bond (HB, SB, or HP) matrices from Simulaid log
                call compare_bondmat(resnames,nrescol,itemp1,temp,nres,
                call askyn(
    c           Compare two residue RMSF lists from Simulaid log
                call askyn(
                call compare_rmsf(resnames,nrescol,siglev,temp,analfile4,
                    cv(ia)=temp(ir)
                call openfile(iw4,0,'Avg difference labeled PDB',26,
                call writeconf(iw4,inpcrdtyp,iobpdb,inpcrdtyporg,
                close (iw4)
    c         Set up adjacency matrix-based analysis
              call quiz(analtyp,idtyp,' ',' ',0,
                call modrepats
                call askyn(
              call getreal('Threshold distance',18,resdistdef,resdistlim,1,
              call quiz(analtyp,iadjtyp,' ',' ',0,'adjacency analysis',18,
              call getint('Highest exponent to raise the adjacency matrix',
              call getint('Power interval to plot',22,1,1,100,npint,0)
              call askyn(
              call getresrange(nsegm,indexs,isegno,ixres,iresno,ifres,
              call askyn('Do you want to mark residues',28,1,0,imarkres,109,
                call openfile(iw4,0,'residue mark',12,'old',analfile4,
              call openps(iw1,xm,ym,title,ltitle76,
    c         Set up distance measuring
              call getint('First  atom number',18,1,1,n,ia1,0)
              call getint('Second atom number',18,1,1,n,ia2,0)
                  call askyn(
                call readax('Axis normal to the plane (1,2,3)',32,3,idax,
                call trnsfr(caref,cres(1,ia1),3)
              call setpbccell('Do you want to use periodic images',34,
                  call pbcdist(caref,c(1,ia2),ia1,ia2,cell,ncell,-iw0,
                  call zeroit(drimg,3)
                  call arrdiff(drimg,cell(1,img),drimg,3)
    c         Set up check for unphysical features
              call getreal('CTFAC',5,1.4,ctfac,1,62)
              call getreal('MINFAC',6,0.4,bondminfac,1,62)
              call getint('MAXDIST',7,min0(nslt,50),1,nslt,maxdist,62)
                call askyn('Do you have LES segments',24,1,-1,iles,0,0)
                    call askyn(qq,27,1,-1,idelseg(is),0,0)
              call setmolres(ifres,ilres,isegno,molresflag,
              call setpbccell('Do you want to use periodic images',34,
              call printbondthres(ialist,nanos,ctfac,bondminfac,
                call prtcell(ioppbc,edge,edge_gen,r,vol,nw,-iw0)
    c         Psudorot calc
                call getring(line,index,ix5,irescol1,irescol2,inamcol1,
                call findrange(isegno,1,nslt,isegpkhx,ifseghx,ilseghx,
                call getint('Residue number of the kink (Proline)',36,
                call findresnum(iresno,ixres,irespro,ifseghx,ilseghx,ia,
              call leftadjustn(resnam,rn,8)
                call askyn('Do you want to project the proline to a plane',
              call getint('Number of helix residues after  the kink',40,
              call findprotbackbone(line,index,iresno,ia,inamcol1,
    c         icab(1)=icapr
    c         icb(1)=icpr
    c         inb(1)=inpr
                call findprotbackbone(line,index,iresno,ia,inamcol1,
                call findprotbackbone(line,index,iresno,ia,inamcol1,
                call getring(line,index,ix5,irescol1,irescol2,inamcol1,
    c           Open Prokink dial window, initialize trajectory accumulators
    c           Helix axis directions
                  call findrange(isegno,1,nslt,iseghx(ihx),ifseghx,ilseghx,
                  call findresnum(iresno,ixres,ir0,ifseghx,ilseghx,ia,
                  call findprotbackbone(line,index,iresno,ia,
              call askyn('Is this a transmembrane protein',31,1,-1,itmem,0,
    c         call getint('First helix of helix pair to track',34,1,1,nhx,
    c    -      ihxt1,0)
    c         call getint('Secnd helix of helix pair to track',34,ihxt1+1,1,
    c    -      nhx,ihxt2,0)
              call askyn('Do you want debug output',24,1,-1,idebughx,0,0)
                call askyn('Do you want to skip bend analysis',33,1,-1,ibx,
                call askyn(linein,65,1,1,isubcrm,0,0)
                call askyn(linein,60,1,1,ioverlay,0,41)
                  call askyn(linein,54,1,1,ioverlaym,0,0)
                  call askyn(linein,59,1,1,ireorienthx,0,0)
                call zeroiti(nhelixok,0,3)
    c         Set up hydropathy labeling
    c         Set up Delphi potential labeling
              call openfile(50,0,'Delphi potential map file',25,'old',
              call readmap(50,xstart,ystart,zstart,gx,gy,gz,ngx,ngy,ngz,
                call getint('Increment to use in reading the potential map',
              call askyn(
                call getreal('Minimum distance between an atom and a grid',
                call askyn(
              call askyn('Do you want to query the potential map',38,1,-1,
    c             Open 2nd grid PDB file
                  call openfile(iw2,0,'grid',4,'new',analfile2,namleno1,
                    call openfile(iw3,0,'grid',4,'new',analfile3,namleno1,
                call delphigrid(iw1,iw2,iw3,c,n,nslt,xstart,ystart,zstart,
    c         Set up circular variance labeling calculation
              call getrange(nsltref_f,1,nsltref_l,nslt,incr,0,
              call getreal(
                call askyn('Do you want to sort solvents by CV',34,1,+1,
                  call getreal('CV threshold to count the # of solvents',
    c         Set up circular variance map calculation
              call openps(iw1,xm,ym-20,title,ltitle76,
    c         Set up DSSP calculation
                call getrange(ifrdssp,1,ilrdssp,nresslt,incr,0,
                  call askyn(
                call askyn(
                call trnsfi(indexdel,ixresno,nresslt)
              call zeroiti(idistdssp,0,9*maxrsd)
    c         Set up RMSD calculations
              call zeroiti(indexdel,0,n)
              call indexit(indexov,1,nslt,0)
                  call askyn(
                  call getreal(
                  call askyn('Do you want to select atoms for overlay',39,
                    call select(line,nrecdel,idcol,asterisk,n,nslt,index,
    c               indexdel: 0 for atoms to use for the overlay, 1 for the rest
                    call masktolist(indexov,indexdel,nslt,nfinalov,0)
    c           indexov contains the list of atoms used for overlay
                call quiz(ans,icalctyp,'a',' ',0,
                  call indexit(indexrmsd,1,nslt,0)
                  call trnsfi(indexrmsd,indexov,nfinalov)
                  call zeroiti(indexdel,0,n)
                  call select(line,nrecdel,idcol,asterisk,n,nslt,index,
                  call masktolist(indexrmsd,indexdel,nslt,nfinalrmsd,0)
    c           indexrmsd contains the list of atoms used for overlay
                call askyn('Do you want to use the input structure instead',
                call askyn('Do you want to plot also the maximum deviation',
                call askyn('Do you want to plot without overlay',35,1,-1,
                call getreal(
                call getreal(
                  call zeroitd(rmsfsum,nresslt)
                  call getreal(
                  call getint('Frequency of plotting error bars',32,
                call zeroitd(cdp,3*nfinalrmsd)
                call zeroitd(cdp2,3*nfinalrmsd)
    c           2D RMSD
                call askyn(
                  call askyn('Do you want RMSD-based clustering',33,
                  call askyn('Do you want to skip plotting the matrix',39,
    c             Read 2-d RMSD map, run clustering
                  call askyn(question,lq,1,1,ians,00,0)
                  call openfile(iw4,0,'previously written .rd2',23,
                  call read_2drmsd(iw4,system,lsystem,trajnam,ltrajnam,
                  close(iw4)
    c             Run clustering
                  call indexit(iconfsel,1,MAX2D,0)
                  call indexit(ixclst,1,MAX2D,0)
                  call askyn(
                  call clusterdistr(nframe,iw0,rmsdlim,rmsdmn,rmsdmx,
                  call trnsfi(ixshuffle,iconfsel,MAX2D)
    c             iconfsel contains the sorted cluster members
    c             indexn, indexo contain the cluster limits
                  call countsim(indexn,indexo,iconfsel,ncl,rdclust,rmsdsim,
                  call askyn('Do you want to replot the input RMSD map',40,
                    call getint(askcolcode,24,maxcolcode,1,maxcolcode,
    c               Generate plot file name
                    call openfile(iw1,0,'(shuffled) RMSD map',19,
                    call openps(iw1,xm_2d,ym_2d,title(1:76),76,pstitle,
                      call indexit(it3,1,MAX2d,0)
                      call adjust_xtraj(xtraj,ifirsttraj,ilasttraj,
                      call plot2drmsd(nrep,iw0,iw1,xtraj,maxrec,title,
                    call indexit(ixshuffle,1,MAX2D,0)
                      call indexit(it3,1,MAX2d,0)
                      call plothead(iw1,xm_2d,ym_2d,title(1:76),76,
                      call plot2drmsd(nrep,iw0,iw1,xtraj,maxrec,title,
                call indexit(iconfsel,1,MAX2D,0)
    c           Cross RMSD
                call getreal(
                  call askyn('Do you want to find matching structures',39,
                  call getreal('MAXimum RMSD for similarity statistics',38,
                  call indexit(iconfsel,1,MAX2D,0)
    c             Read RMSD maps and try to match them
                    call getname(analfile4,namlen4,
                    call openfile(iw4,0,'previously written .rd2',23,
                  call read_2drmsd(iw4,system,lsystem,trajnam1,ltrajnam1,
                  close (iw4)
    c             Run clustering
                  call indexit(index2d1,1,MAX2D,0)
                  call indexit(ixclst,1,MAX2D,0)
                  call clusterdistr(nframe1,iw0,rmsdlim,rmsdmn,rmsdmx,
                  call trnsfi(irepmx1,irepmx,MAX2D)
                  call trnsfi(ixshuffle,index2d1,MAX2D)
                  call countsim(ifclst1,ilclst1,index2d1,ncl1,rdclust,
                    call getname(analfile4,namlen4,
                    call openfile(iw4,0,'previously written .rd2',23,
                  call read_2drmsd(iw4,system,lsystem,trajnam2,ltrajnam2,
                  close (iw4)
    c             Run clustering
                  call indexit(index2d2,1,MAX2D,0)
                  call indexit(ixclst,1,MAX2D,0)
                  call clusterdistr(nframe2,iw0,rmsdlim,rmsdmn,rmsdmx,
                  call trnsfi(irepmx2,irepmx,MAX2D)
                  call trnsfi(ixshuffleref,index2d2,MAX2D)
                  call countsim(ifclst2,ilclst2,index2d2,ncl2,rmsdsim1,
                    call getname(analfile4,namlen4,
                    call openfile(iw4,0,'previously written .rdx',23,
                  call read_2drmsd(iw4,system,lsystem,trajnamr1,ltrajnamr1,
                  call askyn(
    c               Generate plot file name
                    call getint(askcolcode,24,maxcolcode,1,maxcolcode,
                    call openfile(iw1,0,'shuffled cross-RMSD map',23,
                    call openps(iw1,xm_2d,ym_2d,title(1:76),76,pstitle,
                    call adjust_xtraj(xtraj,ifirst,ilast,increment,
                    call plot2drmsd(nrep,iw0,iw1,xtraj,maxrec,title,
                    call indexit(ixshuffle,1,MAX2D,0)
                  call countsimx(ifclst1,ilclst1,index2d1,ncl1,
                  call getreal('Maximum RMSD for mapping',24,rmsdsim1,
                  call mapclustx(2,ifclst1,ilclst1,irepmx1,ncl1,
                  call mapclustx(1,ifclst2,ilclst2,irepmx2,ncl2,
    c         Residue correlation matrix calculation
              call askyn(
              call askyn('Do you also want to plot the covariance matrix',
              call modrepats
                call findat(indexa(ncorr),ifres(ir),ilres(ir),line,index,
    c         Input covariance matrix to get eigenvectors, eigenvalues
              call askyn('Is the covariance matrix file binary',36,1,-1,
                call askyn('Is the matrix broken into 5 numbers/line',40,
                call openfile(iw1,0,'input',5,'old',analfile1,lanalfile1,
              call openfile(iw0,0,'output',6,'new',analfile,lanalfile,
              call askyn('Do you want annotated output',28,1,1,iannout,128,
              call normalmodes(ncorr,iw1,0,iformcov,iw0,
    c         Atom-atom distance distribution calculation
              call askyn(
                call getlist(listpairdist,npairs,1,nslt,2,MAXDDISTR)
                call getclusterpairs(npairs,iclustermem,ifstclst1,ifstclst2,
              call zeroit(pairdistsum,npairs)
              call zeroit(pairdistwsum,2*npairs)
              call zeroit(pairdistsum2,npairs)
              call zeroiti(npairdist,0,10*npairs)
              call getreal('Maximum distance to calculate distribution',42,
    c         Atom-atom distance SD calculation
              call select(line,nrecdel,idcol,asterisk,n,nslt,index,
    c           indexdel: 0 for atoms to use, 1 for the rest
              call masktolist(ianchor,indexdel,nslt,nanchor,0)
              call quiz(ans,isdtyp,'n',' ',0,'SD normalization',16,0,5,6,0)
              call openps(iw1,xm_2d,ym_2d,title(1:76),76,pstitle,
    c         Solvation shell volume calculation
              call getreal('Solvent radius',14,1.4,rsolv,1,0)
              call getint('Number of random points for volume calculation',
              call getint('Print level',11,0,0,0,levout,0)
    c         Principal axis calculation
              call select(line,nrecdel,idcol,asterisk,n,nslt,index,ixres,
              call masktolist(indexa,indexdel,n,natspax,0)
              call condenselist(indexa,natspax,ic,6)
              call condenselist(indexa,natspax,ic,iw0)
                call askyn('Do you want to use the input structure instead',
    c           pstitle='Principal axis evolution'
    c           lpstitle=24
    c         Radius of gyration, hydrodynamic radius, com, dipole moment
              call select(line,nrecdel,idcol,asterisk,n,nslt,index,ixres,
              call masktolist(indexa,indexdel,n,ngyrats,0)
              call condenselist(indexa,ngyrats,ic,6)
    c           pstitle='Gyration and hydrodynamic radii evolution'
    c           lpstitle=41
              call condenselist(indexa,ngyrats,ic,iw0)
    c         Filter solvents
              call getint('Representative atom (center) of the solvent',
              call getreal('Solvent space grid spacing',26,8.0,spacing,1,
              cvmin=0.0
              call quiz(ans,ifilttyp,'c',' ',0,'solvent filtering option',
                  call quiz(ans,intftyp,' ',' ',0,'interface partners',18,
                    call getint('First molecule of the interface',31,1,1,
                    call getint('Second molecule of the interface',32,
                    call getint('Selected interface molecule',27,1,1,
    c               All pairs
                cvrlim=3.5
                cvmin=0.60
                call askyn('Do you want to change these limits',34,1,-1,
                  call getreal('Minimum CV of an interface solvent',34,
                  call getreal(
                  call getreal(
                  call getreal(
                  call getreal(
                  call getreal('Radius of sphere that can not be empty',
    c           Hydrogen bond bridging solvent
                call set_hbondlims(hbf0,hblimfac,angm0,angmin,iw0)
              call strip_cext(analfile1,namleni,namleno,lenext)
              crdexti(1:lenext)=analfile1(namleni-lenext+1:namleni)
    c         print *,'NAMLENI,NAMLENO,LENEXT=',namleni,namleno,lenext
                call askyn(
    c         pstitle='Principal axis calculation'
    c         lpstitle=26
    c!!       Early openps and landscape conflict!
    c       nconfig > 1
          call cofms(cres,crmslt0,nslt,atw)
            call askyn('Do you want mass-weighting',26,1,1,masswt,0,0)
    c       Analysis of single conformations
    c---------(S) Neighbor, bond, angle and torsion list
              call printbondlist(iangpr,itorpr,1,nslt,nslt,c,nneig,ineig,
    c---------(S) 1-4 statistics
              call stat14(c,n,nneig,ineig,index,line,nconfig,pi,iatnm2,
              call findfg(1,n,iatnum,mmtype,indexo,indexn,nhbneig,ih,ifgtyp,
              call findbackbone(1,nslt,iw1,nneig,ineig,nneiga,
    c---------(S) Bond length statistics
              call bondlenstat(c,n,iatnum,nneig,ineig,iatnm2,nconfig,
    c---------(S) Hydrogen-bond list
              call hblist(c,n,nslt,islvw,naslv,iatnum,ifchrg,isegno,nneig,
    c---------(S) Hydrophobic bond or heavy atom VdW contact list
              call nnlisthph_sltb(nslt,ianchor2,iselfanc,indexa,
              call hph_sltblist(ibondtype,c,nslt,nhbneig,ineig,iw0,is1,
    c---------(S) Mutually proximal heavy atom contact list
              call nnlistmpx(n,nanchorr,nanchorn,indexa,indexov,iatnum,
              call mpxblist(n,itemp1,itemp2,cv,itemp3,itemp4,is1,is2,
    c---------(S) Salt bridge list
              call nnlisthph_sltb(nslt,ianchor2,iselfanc,indexa,
              call hph_sltblist(ibondtype,c,nslt,nhbneig,ineig,iw0,
    c---------(S) Calculate residue-residue distances
              call rrdist(c,n,nres,ixresno,iresno,ixres,resnames,ifres,
    c---------(S) Calculate a distance
              call pbcdist(c(1,ia1),c(1,ia2),ia1,ia2,cell,ncell,0,0,
    c---------(S) Check for potentially unphysical contacts
              call checkunphys(c,nslt,n,naslv,islvw,iatnum,ifchrg,isegno,
    c---------(S) Pseudorotation angle calculation
                call getring(line,index,ix5,irescol1,irescol2,inamcol1,
                call pseudorot(c,n,ix5,nmem,nneig,ineig,psr5,psr,npsr,
    c---------(S) Calculate Proline kinks
              call prokinkcalcla(1,c,nslt,bend,wobble,faceshift,
    c---------(S) Hydropathy labeling
              call hydropathylist(n,nslt,ixres,resnames,cv,ihydtyp,
              call writeconf(iw0,inpcrdtyp,inpcrdtyp,inpcrdtyporg,n,n,nslt,
    c---------(S) Circular variance labeling of solute atoms and solvent molecules
              call cvlist(c,n,nslt,nsltref_f,nsltref_l,naslv,islvrep,icvtyp,
              call writeconf(iw0,inpcrdtyp,inpcrdtyp,inpcrdtyporg,n,n,nslt,
    c---------(S) Circular variance residue-residue plot
              call cvplot(c,n,nslt,icvtyp,line,index,indexs,inamcol1,
    c---------(S) DSSP secondary structure assignment
              call dssp(c,1,n,nslt,line,index,inamcol1,inamcol2,iresncol1,
    c---------(S) Hydrogen-bond bridge analysis
              call hblist(c,n,nslt,islvw,naslv,iatnum,ifchrg,isegno,nneig,
              call hbbridge(nanchor,ianchor,indexa,ianchor2,iselfanc,lpath,
              call hbbridgeprint(nanchor,ianchor,lpath,nbridgetype,
    c---------(S) Ramachandran plot
              call ramachandran(c,nslt,index,line,nconfig,pi,nresfound,
              call ramachandranplot(nresfound,iw1,xm,iallrama)
    c---------(S) Angle dial plots
              call angledials(c,nslt,nangsel,ixtor1234,index,line,pi,
    c---------(S) Torsion dial plots
              call torsiondials(c,nslt,ntorsel,ixtor1234,index,line,pi,
    c---------(S) Delphi map annotation
              call delphilabel(c,n,nslt,xstart,ystart,zstart,gx,gy,gz,cv)
              call writeconf(iw0,inpcrdtyp,inpcrdtyp,inpcrdtyporg,n,n,nslt,
    c---------(S) Helix axis directions
                call helixaxis(c,nslt,iw0,calph(1,1,ihx),axisdir(1,ihx),
                call dssp(c,iadssp1,n,iadssp2,line,index,inamcol1,
                call checkforhelix(ireshx2(ihx)-ireshx1(ihx)+1,
    c---------(S) Distance calculation
                call pairdistcalc(c,nslt,npairs,listpairdist,pairdistsum,
                call clusterdistcalc(c,nslt,npairs,iclustermem,
              call pairdistprint(nframe,npairs,listpairdist,iclusterdist,
    c---------(S) Solvation shell volume calculation
              call volcalc(nrand,c,nslt,isegno,iatnum,nsegslt,molsltlim,
    c---------(S) Principal axis calculation
              call princax(c,c2,atw,temp,n,indexdel,evecs0,evals0,iw0,0,1,
              call askyn('Do you want to align the structure to an axis',
                call getint('Axis index (1/2/3) to align',27,3,1,3,iax,000)
                call cofms(c,crmslt,nslt,aw)
                call shiftmol(c,nslt,crmslt,c2,-1.0)
                call rotate_c(c2,n,rot,c2,'PRINCAX',7)
                call writeconf(iw1,inpcrdtyp,inpcrdtyp,inpcrdtyporg,n,n,
                close (iw1,status='delete')
    c---------(S) Radius and dipole calculation
                call molrad(c,indexa,ngyrats,iw0,MAXREC)
                call celldipole(c,n,nslt,indexa,ngyrats,charge,icharges,atw,
    c---------(S) Calculate adjacency-matrix based analysis
              call rrconn(c,n,ires1,ires2,ifres,ilres,iatnum,ignoreh,
    c---------(S) Calculate molecule-molecule distance matrix
                call mmdist(c,n,atw,iatnum,nmolslt,molsltlim,c2,temp,
    c---------(S) Filter solvents
              call filterslv(c,iatnum,nslt,n,naslv,numsolv,molsltlim,
                  call writeint(analfile1,namleno,nconfig,lenc)
                  call openfile(45,0,'FILT',4,'new',analfile1,namleno,
                call writeout(45,inpcrdtyp,inpcrdtyp,line,index,isegno,
                close (45)
    c       Trajectory scan
    c       Open, initialize trajectory
    c       2D RMSD calculation in blocks
              call opentraj(c,1,inpt,inptrajtyp,n,ntitltr,trtitle,
            call save_traj_lim(ifirst,ilast,increment,1)
            call adjust_xtraj(xtraj,ifirst,ilast,increment,framefac,
    c         Read a conformation
              call readtraj(inpt,inptrajtyp,mmctrajtyp,n,naslv,nwatr,
                call selectconf(numsel,ninconf,ifirst,increment,
    c           Calculate RMSD between the ncop structures in c
                    call rmsd(c(1,ic2),c(1,ic1),nslt,nfinalov,nfinalrmsd,
                    call progress_rep(0,nframe2d,nframesign)
    c           print *,'Self block done nframe1=',nframe1
    c           Open same trajectory and skip to nframeread+increment
                call opentraj(c2,0,inpt,inptrajtyp,n,ntitltr,trtitle,
    c           Read blocks into cref and calculate RMSD with the c block
    c             Read a conformation
                  call readtraj(inpt,inptrajtyp,mmctrajtyp,n,naslv,nwatr,
                    call selectconf(numsel,ninconf2,ifirst,increment,
    c               RMSD calc between c and cres
                        call rmsd(c(1,ic1),cres(1,ic2),nslt,nfinalov,
                        call progress_rep(0,nframe2d,nframesign)
    c               print *,'Block pair done nframe1,2=',nframe1,nframe2
                  call opentraj(c,1,inpt,inptrajtyp,n,ntitltr,trtitle,
                    call readtraj(inpt,inptrajtyp,mmctrajtyp,n,naslv,nwatr,
            close (inpt)
              call openps(iw1,xm,ym,title(1:76),76,pstitle,lpstitle,trajnam,
    c       Cross-RMSD in blocks
              call opentraj(c,1,inpt,inptrajtyp,n,ntitltr,trtitle,
            call save_traj_lim(ifirst1,ilast1,increment1,1)
              call askstop(0)
            call adjust_xtraj(xtraj,ifirst1,ilast1,increment1,framefac,
    c         Read a conformation
              call readtraj(inpt,inptrajtyp,mmctrajtyp,n,naslv,nwatr,
                call selectconf(numsel,ninconf1,ifirst1,increment1,
    c           Read blocks of traj2
                call opentraj(cres,0,inpt2,inptrajtyp,n,ntitltr,trtitle,
                call save_traj_lim(ifirst2,ilast2,increment2,2)
                  call askstop(0)
    c             Read a conformation
                  call readtraj(inpt2,inptrajtyp,mmctrajtyp,n,naslv,nwatr,
                    call selectconf(0,ninconf2,ifirst2,increment2,
    c               RMSD calc between c and cres
    c               print *,'NFRAME1,NCOP1,NFRAME2,NCOP2=',
    c    -            NFRAME1,NCOP1,NFRAME2,NCOP2
                        call rmsd(c(1,ic1),cres(1,ic2),nslt,nfinalov,
                        call progress_rep(0,nframe2d,nframesign)
    c               print *,'NINCONF1,2=',ninconf1,ninconf2,
    c    -            ' ILAST1,2=',ilast1,ilast2
                close (inpt2)
              call openps(iw1,xm,ym,title(1:76),76,pstitle,lpstitle,trajnam,
            call opentraj(c,nrep+nrepinc,inpt,inptrajtyp,n,ntitltr,trtitle,
            call save_traj_lim(ifirst,ilast,increment,1)
            call adjust_xtraj(xtraj,ifirst,ilast,increment,framefac,
    c         Round of the max # of frames for plots
              call roundlimint(maxtrajplot,idiv,ndivdssp)
              call lastchar(title,ltitle,80)
    c         PS files with trajectory name
                call getint('Last configuration in the trajectory',36,
                  close (inpt)
                call askstop(0)
                    call helixaxis(cres,nslt,iw0,calph0(1,1,ihx),
    c               Establish helix direction(s)
    c                 Establish intra/extra cellular helix end labeling
                      call quiz(ans,ihxtyp,tmchar,' ',0,'Helix type',10,0,
    c               Calculate reference turn angles
                      call angcomp(perpvec0(1,ir-1,ihx),axisdir0(1,ihx),
    c               Shift the helix so that the start is at the origin
                      call dvdif(calph0(1,ir,ihx),axisini0(1,ihx),
    c       Open, read and analyze trajectory
    cx      print *,'Start scan nmc=',nmc
    c         Read a conformation
              call readtraj(inpt,inptrajtyp,mmctrajtyp,n,naslv,nwatr,
                call askstop(-1)
    c         if (ianaltyp .eq. 23) irepscan=0
                call selectconf(numsel,ninconf,ifirst,increment,
    c           print *,'Selected nmc=',nmc
                      call pbcdist(c(1,ia1),c(1,ia2),ia1,ia2,cell,ncell,
                      call zeroit(drimg,3)
                      call arrdiff(drimg,cell(1,img),drimg,3)
    c               Update box info (assuming isotropic box fluctuation)
                    call updatecell(inptrajtyp,edge)
                  call checknnlist(1,n,ineig,nneig,nerr,maxng)
                  call comparetop(c,n,nneig,ineig,iatnum,innlist,nslt,naslv,
    c           Configuration is ready, choose the analysis
    c-------------(T) Neighbor, bond, angle and torsion list
                  call printbondlist(iangpr,itorpr,1,nslt,nslt,c,nneig,
    c-------------(T) 1-4 statistics
                  call stat14(c,n,nneig,ineig,index,line,nframe,pi,iatnm2,
                  call findfg(1,n,iatnum,mmtype,indexo,indexn,nhbneig,ih,
    c-------------(T) Bond length statistics
                  call bondlenstat(c,n,iatnum,nneig,ineig,iatnm2,nframe,
    c-------------(T) Hydrogen-bond list
                  call hblist(c,n,nslt,islvw,naslv,iatnum,ifchrg,isegno,
                    call selectbond(ixres,nanchor,ianchor,indexa,ianchor2,
    c-------------(T) Hydrophobic bond or heavy atom contact list
                  call nnlisthph_sltb(nslt,ianchor2,iselfanc,indexa,
                    call selectbond(ixres,nanchor,ianchor,indexa,ianchor2,
    c-----------(T) Mutually proximal heavy atom contact list
                call nnlistmpx(n,nanchorr,nanchorn,indexa,indexov,iatnum,
                call mpxblist(n,itemp1,itemp2,cv,itemp3,itemp4,is1,is2,
    c-------------(T) Salt bridge list
                  call nnlisthph_sltb(nslt,ianchor2,iselfanc,indexa,
                    call selectbond(ixres,nanchor,ianchor,indexa,ianchor2,
    c-------------(T) Calculate residue-residue distances
                  call rrdist(c,n,nres,ixresno,iresno,ixres,resnames,ifres,
    c-------------(T) Calculate a distance
    c                 See if c(1,ia1) has not been reset to the cell
                      call pbcdist(caprev1,c(1,ia1),ia1,ia2,cell,ncell,-iw0,
    c                   PBC change found
                        call arrsum(drimg,cell(1,img),drimg,3)
    c                 See if c(1,ia2) has not been reset to the cell
                      call pbcdist(caprev2,c(1,ia2),ia1,ia2,cell,ncell,-iw0,
    c                   PBC change found
                        call arrdiff(drimg,cell(1,img),drimg,3)
                    call trnsfr(caprev1,c(1,ia1),3)
                    call trnsfr(caprev2,c(1,ia2),3)
                    call arrdiff(c(1,ia2),c(1,ia1),dr,3)
                    call arrsum(dr,drimg,dr,3)
    c                 See if c(1,ia2) has not been reset to the cell
                      call pbcdist(caprev2,c(1,ia2),ia1,ia2,cell,ncell,-iw0,
    c                   PBC change found
                        call arrdiff(drimg,cell(1,img),drimg,3)
                    call trnsfr(caprev2,c(1,ia2),3)
                    call arrdiff(c(1,ia2),caref,dr,3)
                    call arrsum(dr,drimg,dr,3)
                  call trajlimtest(nframe,MAXFRAMES)
    c-------------(T) Check for potentially unphysical contacts
                  call checkunphys(c,nslt,n,naslv,islvw,iatnum,ifchrg,
    c-------------(T) Pseudorotation angle calculation
                  call pseudorot(c,n,ix5,nmem,nneig,ineig,psr5,psr,npsr,
    c-------------(T) Calculate Proline kinks
                  call prokinkcalcla(nrep,c,nslt,bend,wobble,
                  call pseudorot(c,n,ix5,nmem,nneig,ineig,psr5,psr,npsr,
    c-------------(T) Hydropathy labeling
                  call hydropathylist(n,nslt,ixres,resnames,cv,ihydtyp,
                  call writeconf(iw0,inpcrdtyp,inpcrdtyp,inpcrdtyporg,n,n,
    c-------------(T) Circular variance labeling of solute ats and solvent molecs
                  call cvlist(c,n,nslt,nsltref_f,nsltref_l,naslv,islvrep,
                  call writeconf(iw0,inpcrdtyp,inpcrdtyp,inpcrdtyporg,n,n,
    c-------------(T) DSSP secondary structure assignment
                  call dssp(c,1,n,nslt,line,index,inamcol1,inamcol2,
                  call plotdssp(iw1,ifsse,ilsse,itypsse,nsse,ninconf,
    c-------------(T) Hydrogen-bond bridge analysis
                  call hblist(c,n,nslt,islvw,naslv,iatnum,ifchrg,isegno,
                  call hbbridge(nanchor,ianchor,indexa,ianchor2,iselfanc,
    c-------------(T) Ramachandran plot
                  call ramachandran(c,nslt,index,line,nframe,pi,nresfound,
                  call ramachandranplot(nresfound,iw1,xm,iallrama)
    c-------------(T) Torsion dial plots
                call torsiondials(c,nslt,ntorsel,ixtor1234,index,line,pi,
    c-------------(T) Delphi map labeling
                  call delphilabel(c,n,nslt,xstart,ystart,zstart,gx,gy,gz,
                  call writeconf(iw0,inpcrdtyp,inpcrdtyp,inpcrdtyporg,n,n,
    c-------------(T) Helix directions
    cd77              write (77,*) ' ------------ Nframe=',nframe,' -----------'
                      call dssp(cres,iadssp1,n,iadssp2,line,index,
                      call checkforhelix(ireshx2(ihx)-ireshx1(ihx)+1,
                    call helixcomp(nslt,nreshx(ihx),calph0(1,1,ihx),
                      call dssp(c,iadssp1,n,iadssp2,line,index,inamcol1,
                      call checkforhelix(ireshx2(ihx)-ireshx1(ihx)+1,
    c-------------(T) 1-D RMSD and residue RMS fluctuations
                  call rmsd(cres,c,nslt,nfinalov,nfinalrmsd,atw,atwsum,temp,
                      cdp(k,ia)=cdp(k,ia)+c2(k,indexrmsd(ia))
                      cdp2(k,ia)=cdp2(k,ia)+c2(k,indexrmsd(ia))**2
    c               write (iw0,*) ia,' indexrmsd(ia)=',indexrmsd(ia)
                    call rmsf(chn,c2,atw,nfinalrmsd,ixres,indexrmsd,
    c               write (77,8791) nframe,(rmsfsum(i),i=1,nresslt)
    c8791           format(' NFRAME=',i4,' rmsfsum:',/,(10f8.2))
    c-------------2D RMSD and cross RMSD calculations done separately (in blocks)
    c-------------(T) Residue correlation calculation
                  call residcorr(c,c,nslt,indexa,ncorr,nframe)
    c-------------(T) Atom-atom distribution calculation
                    call pairdistcalc(c,nslt,npairs,listpairdist,
                    call clusterdistcalc(c,nslt,npairs,iclustermem,
    c-------------(T) Solvation shell calculation
                  call volcalc(nrand,c,nslt,isegno,iatnum,nsegslt,molsltlim,
                  call trajlimtest(nframe,MAXFRAMES)
    c-------------(T) Principal axis calculation
                  call princax(c,c2,atw,temp,n,indexdel,evecs0,evals0,iw0,
    c-------------(T) Radius and dipole calculation
                  call molrad(c,indexa,ngyrats,iw0,MAXREC)
                  call celldipole(c,n,nslt,indexa,ngyrats,charge,icharges,
    c-------------(T) Angle dial plots
                  call angledials(c,nslt,nangsel,ixtor1234,index,line,pi,
    c-------------(T) Calculate molecule-molecule distance matrix
                    call mmdist(c,n,atw,iatnum,nmolslt,molsltlim,c2,temp,
    c-------------(T) Atom-atom distance SD calculation
                  call atomdist_sd(c,nslt,ianchor,nanchor,nframe)
    c-------------(T) Solvent filtering
                  call filterslv(c,iatnum,nslt,n,naslv,numsolv,molsltlim,
                      call writeint(analfile,namleno,nframe,lenc)
                      call openfile(45,0,'FILT',4,'new',analfile1,namleno,
                    call writeout(45,inpcrdtyp,inpcrdtyp,line,index,isegno,
    c             Print progress report
                    call progress_rep(nframe,nframe2d,nframesign)
    c               if (ianaltyp .eq. 23)
    c    -            nframetot=nframetot*(nframetot-1)/2
    c       Trajectory scan done
            call testconst(0,1,2,0.0,1.0,2.0,6,nfail,1,'ANPO')
    c       Post scan activities
    c         Use trajectory title
              call dialps(iw1,prokinklab,lprokinklab,title,ltitle,resrange,
              call trajstat(iw0,ndials,5,prokinklab,lprokinklab,0,1,0,
                call arminmax2(res(1,1,1),1,nframe,2,armin1,armax1,
                call roundlim(armax1,y1div,ny1div)
                call plot2fun(iw1,2,xtraj,res(1,1,1),res(1,1,1),nframe,
                call arminmax2(res(1,1,3),1,nframe,1,armin1,armax1,
                call roundlim(armax1,y1div,ny1div)
                call plot2fun(iw1,1,xtraj,res(1,1,3),res(1,1,3),nframe,
                call arminmax2(res(1,1,2),1,nframe,2,armin1,armax1,
                call roundlim(armax1-xmn,xdv,nxdv)
                call roundlim(armax2-ymn,ydv,nydv)
                call plot2fun(iw1,2,xtraj,res(1,1,2),res(1,1,2),nframe,
                call plot2d(iw1,res(1,1,2),nframe,nfravgt,xmn,xdv,nxdv,
                  call dialps(iw1,helixang,lhelixang,title,ltitle,resrange,
                  call trajstat(iw0,ndials,6,helixang,lhelixang,10,10,7,
                  call arminmax2(res(1,1,incrhx+7),1,nframe,2,armin1,armax1,
                  call roundlim(armax2-y2mn,y2div,ny2div)
                  call plot2fun(iw1,2,xtraj,res(1,1,incrhx+7),
                  call arminmax2(res(1,1,incrhx+8),1,nframe,2,armin1,armax1,
                  call roundlim(armax1-y1mn,y1div,ny1div)
                  call roundlim(armax2-y2mn,y2div,ny2div)
                  call plot2fun(iw1,2,xtraj,res(1,1,incrhx+8),
    c             Plot the progress of the helix center in the plane of the other ax
                  call arminmax2(res(1,1,incrhx+9),1,nframe,2,armin1,armax1,
                  call roundlim(armax1-xmn,xdv,nxdv)
                  call roundlim(armax2-ymn,ydv,nydv)
                  call plot2fun(iw1,2,xtraj,res(1,1,incrhx+9),
                  call plot2d(iw1,res(1,1,incrhx+9),nframe,nfravgt,xmn,
    c             Plot the move of the helix start in the plane of the other axes
                  call arminmax2(res(1,1,incrhx+10),1,nframe,2,armin1,
                  call roundlim(armax1-xmn,xdv,nxdv)
                  call roundlim(armax2-ymn,ydv,nydv)
                  call plot2fun(iw1,2,xtraj,res(1,1,incrhx+10),
                  call plot2d(iw1,res(1,1,incrhx+10),nframe,nfravgt,
    c             Plot the progress of the helix end in the plane of the other axes
                  call arminmax2(res(1,1,incrhx+11),1,nframe,2,armin1,
                  call roundlim(armax1-xmn,xdv,nxdv)
                  call roundlim(armax2-ymn,ydv,nydv)
                  call plot2fun(iw1,2,xtraj,res(1,1,incrhx+11),
                  call plot2d(iw1,res(1,1,incrhx+11),nframe,nfravgt,
    c             Plot the progress of the helix plane normal's projection
                  call arminmax2(res(1,1,incrhx+12),1,nframe,2,armin1,
                  call roundlim(-10.0*xmin,xdv,nxdv)
                  call plot2d(iw1,res(1,1,incrhx+12),nframe,nfravgt,
                    cv(i)=(180./3.141592)*dacoscheck(ddd,scp_ax,0,6,
                  call arminmax2(cv,1,nframe,1,armin1,armax1,armin2,armax2,
                  call roundlim(armax1-y1mn,y1div,ny1div)
                  call plot2fun(iw1,1,xtraj,cv,cv,nframe,0.0,0.0,00,y1mn,
                    call blankout(printrlab(i),1,25)
    c                 Plot the two helix end-end distances
    c                 res(2,nframes,ixres)=10000*idistss+idisthh
                        call trnsfr(res(1,i,1),plotdat(1,i),2)
                      call arminmax2(plotdat,1,nframe,2,armin1,
                      call roundlim(armax2-y2mn,y2div,ny2div)
                      call plot2fun(iw1,2,xtraj,plotdat,plotdat,nframe,0.0,
    c                 Plot the two kinds of helix-helix distances
    c                 res(1,nframes,ixres)=10000*idist+idist_cc
                        call trnsfr(res(1,i,2),plotdat(1,i),2)
                      call arminmax2(plotdat,1,nframe,2,armin1,armax1,
                      call roundlim(armax2-y2mn,y2div,ny2div)
                      call plot2fun(iw1,2,xtraj,plotdat,plotdat,nframe,0.0,
                      call trajstat(iw0,0,1,helixang,lhelixang,4,4,
    c             Plot helix axis angles
    c                 Convert angles back to sin/cos
    c                 res(1,nframes,nhx2tot+ixres)=10000*iang+idhang
                  call dialps(iw1,hxhxlab,lhxhxlab,title,ltitle,resrange,
                  call trajstat(iw0,nhx2tot,nhx2tot,hxhxlab,lhxhxlab,0,1,
    c             Plot helix axis dihedral angles
    c                 Convert angles back to sin/cos
    c                 res(1,nframes,nhx2tot+ixres)=10000*iang+idhang
                  call dialps(iw1,hxhxlab,lhxhxlab,title,ltitle,resrange,
                  call trajstat(iw0,nhx2tot,nhx2tot,hxhxlab,lhxhxlab,0,1,
    c             res(1,nframes,2*nhx2tot+ixres)=10000*idist_ar1+idist_ar2
    c             Plot helix rotations
                  call dialps(iw1,hxhxlab,lhxhxlab,title,ltitle,resrange,
                  call trajstat(iw0,nhx2tot2,nhx2tot2,hxhxlab,lhxhxlab,0,1,
    c           Torsion angle dial plots
                call dialps(iw1,talab,ltalab,title,ltitle,
                call trajstat(iw0,ntorsel,MAXCOPY1,talab,ltalab,0,1,0,
                  call askyn(
                  call askyn(
                    call indexit(index2d,1,ndials,0)
                    call read_write_ccc(ndials,iout_ccc,+1)
                    call normalmodes(ndials,iout_ccc,nframe,1,iw0,1,ierr,
                    call read_write_ccc(ndials,iout_ccc,-1)
                    close (iout_ccc,status='delete')
                  call askyn(
                call arminmax2(res(1,1,7),1,nframe,2,armin1,armax1,
                call roundlim(armax1,y1div,ny1div)
                call arminmax2(res(2,1,7),1,nframe,2,armin1,armax1,
                call roundlim(armax1,y2div,ny2div)
                call plot2fun(iw1,maxdevplot+1,xtraj,res(1,1,7),res(1,1,7),
                  call arminmax2(res(1,1,8),1,nframe,2,armin1,armax1,
                  call roundlim(armax1,y1div,ny1div)
                  call roundlim(armax2,y2div,ny2div)
                  call plot2fun(iw1,maxdevplot+1,xtraj,res(1,1,8),
                call trajstat(iw0,0,6,helixang,lhelixang,4,4,7,rmsdlab,
                  call correl(iw0,res(1,1,7),1,'RMSD',4,av1,sd1,res(1,1,9),
                  call arminmax2(res(1,1,9),1,nframe,2,armin1,armax1,
                  call scatterps(iw1,rmsdmin,rmsdmax,575.0,armin1,armax1,
    c               print *,i,' nframes_err=',nframes_err(i)
                  call rmsf_av(cdp,cdp2,atw,nfinalrmsd,indexrmsd,ixres,
                    call getrange(ifstplot,1,ilstplot,nresslt,incr,0,
                    call blockfromcum(bl,distrerr(1,ir),xcum,MAXDISTRN)
                    call batchmean(MAXDISTRN,0,bl,'RMSF error',10,iw0,1,
                    ci=err12(2,ir)
                  call arminmax2(rmsfavs,ifstplot,ilstplot,2,armin1,
                    call arminmax2(rmsfavs,ifstplot,ilstplot,2,armin1,
                    call roundlim(amax1(armax1,armax2),y1div,ny1div)
                    call askyn('Do you want to use the residue number',37,
                    cv(ir)=ifstres-1+ir
                  call askyn('Do you want error bars plotted on RMSF(Cref)',
                    call roundlim(armax,y1div,ny1div)
                  call plot2fun(iw1,2,cv,rmsfavs(1,ifstplot),
                  call askyn('Do you want to write a PDB file with RMSF',41,
                        cv(ia)=rmsfav(ir)
                    call openfile(iw2,0,'RMSF',4,'new',analfile2,namleno2+4,
                    call writeconf(iw2,inpcrdtyp,inpcrdtyp,inpcrdtyporg,n,n,
    c           Prepare 2-D RMSD plot
                call plot2drmsd(nrep,iw0,iw1,xtraj,maxrec,title,'',0,
                  call clusterdistr(nframe,iw0,rmsdlim,rmsdmn,rmsdmx,
    c             iconfsel contains the sorted cluster members
    c             indexn, indexo contain the cluster limits
                  call countsim(indexn,indexo,iconfsel,ncl,rdclust,rmsdsim,
                    call trnsfi(ixshuffle,iconfsel,MAX2D)
                    call indexit(it3,1,MAX2D,0)
                      call openfile(iw1,0,'analysis',8,'old',analfile1,
                    call plot2drmsd(nrep,iw0,iw1,xtraj,maxrec,title,
                    call askyn('Do you want to inspect the clustered map'//
                      call clusterplot(iw1,xtraj,value,indexn,indexo,ncl,
                      close (iw1)
                  call clusterplot(iw1,xtraj,value,indexn,indexo,ncl,ixclst,
    c           Prepare cross RMSD plot
                call plot2drmsd(nrep,iw0,iw1,xtraj,maxrec,title,'',0,
    c           Prepare distance distribution output
                call pairdistprint(nframe,npairs,listpairdist,iclusterdist,
                call askyn(
                  call openfile(iw4,0,'PDB with distances marked',25,
                  call writeconf(iw4,inpcrdtyp,iobpdb,inpcrdtyporg,
                  close (iw4)
    c           Prepare distance SD matrix plot
                call plot_atomdist_sd(nslt,line,index,inamcol1,inamcol2,
                call trajstat(iw0,0,1,helixang,lhelixang,4,4,1,volumelab,
                call arminmax2(res(1,1,1),1,nframe,2,armin1,armax1,
                call setdivxy(armin1,armax1,ny1div,y1div,y1min)
                call setdivxy(armin2,armax2,ny2div,y2div,y2min)
                call plot2fun(iw1,2,xtraj,res(1,1,1),res(1,1,3),nframe,
                call arminmax2(res(1,1,2),1,nframe,2,armin21,armax21,
                call setdivxy(armin21,armax21,ny1div,y1div,y1min)
                call setdivxy(armin22,armax22,ny2div,y2div,y2min)
                call plot2fun(iw1,2,xtraj,res(1,1,2),res(1,1,4),nframe,
                call zeroiti(indexo,0,nbin)
                call roundlim(range1,xdiv,nxdiv)
                  cv(i)=xmin+(i-0.5)*range1/nbin
                  charge(i)=indexo(i)
                call roundlim(histmax,y1div,ny1div)
                call plot2fun(iw1,1,cv,charge,res(1,1,3),nbin,xmin,xdiv,
                call zeroiti(indexo,0,nbin)
                call setdivxy(armin2,armax2,nxdiv,xdiv,xmin)
                  cv(i)=xmin+(i-0.5)*(armax2-armin2)/nbin
                  charge(i)=indexo(i)
                call roundlim(histmax,y2div,ny2div)
                call plot2fun(iw1,1,cv,charge,res(1,1,3),nbin,xmin,xdiv,
    c           Principal axis calculation
                call arminmax2(res(1,1,6),1,nframe,2,armin1,armax1,
                call setdivxy(armin1,armax1,ny1div,y1div,y1min)
                call setdivxy(armin2,armax2,ny2div,y2div,y2min)
                call setdivxy(0.0,xtraj(nframe),nxdiv,xdiv,xmin)
                call plot2fun(iw1,2,xtraj,res(1,1,6),res(1,1,6),nframe,xmin,
                call arminmax2(res(1,1,7),1,nframe,1,armin1,armax1,
                call setdivxy(armin1,armax1,ny1div,y1div,y1min)
                call plot2fun(iw1,1,xtraj,res(1,1,7),res(1,1,7),nframe,xmin,
    c           Molecular radii, moments of inertia and dipole moment calculation
                call averageres(nframe,res,1,1,MAXFRAMES,MAXCOPY1,rgav,rgsd)
                call averageres(nframe,res,1,2,MAXFRAMES,MAXCOPY1,rhav,rhsd)
                call averageres(nframe,res,2,2,MAXFRAMES,MAXCOPY1,rmxmn,
                call arminmax2(res(1,1,1),1,nframe,2,armin1,armax1,
                call setdivxy(armin1,armax1,ny1div,y1div,y1min)
                call setdivxy(armin2,armax2,ny2div,y2div,y2min)
                call setdivxy(0.0,xtraj(nframe),nxdiv,xdiv,xmin)
                call plot2fun(iw1,2,xtraj,res(1,1,1),res(1,1,1),nframe,xmin,
                  call averageres(nframe,res,1,10,MAXFRAMES,MAXCOPY1,avcell,
                  call averageres(nframe,res,2,10,MAXFRAMES,MAXCOPY1,avslt,
                  call arminmax2(res(1,1,10),1,nframe,2,armin1,armax1,
                  call setdivxy(armin1,armax1,ny1div,y1div,y1min)
                  call setdivxy(armin2,armax2,ny2div,y2div,y2min)
                  call averageres(nframe,res,1,11,MAXFRAMES,MAXCOPY1,
                  call averageres(nframe,res,2,11,MAXFRAMES,MAXCOPY1,
                  call averageres(nframe,res,1,12,MAXFRAMES,MAXCOPY1,
                    call averageres(nframe,res,1,13,MAXFRAMES,MAXCOPY1,
                    call averageres(nframe,res,2,13,MAXFRAMES,MAXCOPY1,
                    call averageres(nframe,res,1,14,MAXFRAMES,MAXCOPY1,
                  call plot2fun(iw1,nplot,xtraj,res(1,1,10),res(1,1,10),
                  call correl(iw0,res(1,1,1),1,'Radius of gyration',18,
    c           Angle dial plots
                call dialps(iw1,talab,ltalab,title,ltitle,
                call trajstat(iw0,nangsel,MAXCOPY1,talab,ltalab,0,1,0,
    c         Print contact count
              call blankout(linein,1,80)
    c         print *,'NSEGCOL=',nsegcol
                  call lastchar(linein,lc,80)
              call print_rrdist(itypavg,nframe,irefres1,irefres2,inegres1,
              call askyn(
                call quiz(ans,icontyp,' ',' ',0,'Contact type',12,0,5,6,
                      cv(ia)=irescount1(ir)
                      cv(ia)=irescount2(ir)
                      cv(ia)=irescount3(ir)
                      cv(ia)=irescount1(ir)+irescount2(ir)+irescount3(ir)
                call openfile(iw4,0,'Contact-labeled PDB',19,
                call writeconf(iw4,inpcrdtyp,iobpdb,inpcrdtyporg,
                close (iw4)
              call dialps(iw1,prokinklab,lprokinklab,title,ltitle,resrange,
              call trajstat(iw0,ndials,5,prokinklab,lprokinklab,0,1,0,
              call openps(iw1,xm,ym-30,title,ltitle,
              call plotresidcorr(ncorr,nframe,indexa,indexs,
    c           Calculate eigenvalues/eigenvectors
                close (iw0)
                call openfile(iw0,0,'eigenvalues & eigenvectors',26,'new',
                call normalmodes(ncorr,iucorrmat,nframe,
                close (iucorrmat,status='delete')
                call rainbowscale(iw1,50,450,25,nframe,xtraj(nframe),0.0,
                call dialps(iw2,ramalab,lramalab,title,ltitle,' ',0,5,
                close (iw2)
                call ramachandran_hist(nresfound,resnames,nrescol,iw0,1,0)
                call ramachandran_hist(nresfound,resnames,nrescol,iw0,0,1)
                call openfile(iw3,0,'2D_Rama_trace',13,'new',analfile3,
                call openfile(iw2,0,'Psi-Phi autoc',13,'new',analfile2,
                  cosphisum=0.d0
                  cospsisum=0.d0
                    cosphisum=cosphisum+res(1,i,2*ir-1)
                    cospsisum=cospsisum+res(1,i,2*ir)
                  cosphisum=cosphisum/nframe
                  cosphisum=cosphisum/sqsum
                  cvphi=1.0-sqsum
                  cospsisum=cospsisum/nframe
                  cospsisum=cospsisum/sqsum
                  cvpsi=1.0-sqsum
                  call plot2d(iw3,xyplot,nframe,1,-180.0,30.0,12,
                  call resautocorr(ir,incrementac,ncf,xyplot,MAXFRAMES)
                  call plot2fun(iw2,1,rplot,xyplot,xyplot,ncf,0.0,0.0,0,
              close (iw1)
    c           Print statistics
              close (iw2)
              close (iw3)
              call hbbridgeprint(nanchor,ianchor,lpath,nbridgetype,
              call finalizebonds(n,nbfound,nbfoundorig,nbresfound,iw0,iw1,
              call finalizebonds(n,nbfound,nbfoundorig,nbresfound,iw0,iw1,
              call finalizebonds(n,nbfound,nbfoundorig,nbresfound,iw0,iw1,
              call finalizebonds(n,nbfound,nbfoundorig,nbresfound,iw0,iw1,
              call finalizebonds(n,nbfound,nbfoundorig,nbresfound,iw0,iw1,
                cs=cospsrs(m)/dfloat(nframe)
                cvm=1.d0-dsqrt(sinpsrs(m)**2+cospsrs(m)**2)/dfloat(nframe)
                call askyn(
                  call getname(trackfile,ltrackfile,
                  call openfile(92,0,' ',1,'new',trackfile,ltrackfile,
                  call writetrack(iout_track,iw0,30,nbfoundorig,nres2d,
              call askyn(question(1:60),60,1,-1,iwpdb,0,0)
                call askyn('Do you want a PDB file with lines between the'//
                call askyn(
                call askyn(
                        cv(ia)=sum
                call openfile(iw4,0,'Percent bond labeled PDB',24,
                call writeconf(iw4,inpcrdtyp,iobpdb,inpcrdtyporg,
                  call writeconn(iw4,ifres,ilres,line,index,inamcol1,
                close (iw4)
          call datprt(iw0,version,0,mark0,lmark0,hostname,lhostname,
          close (iw0)
          call testconst(0,1,2,0.0,1.0,2.0,6,nfail,1,'ANTD')
    c2147  format(' WARNING: x and z cell change factors differ:',2f10.7,/,
    c     -  10x,'Cell change may be anisotropic')
          call getreal('H-bond length tolerance factor',30,hblimfac_def,
          call getreal('A-H...B angle minimum accepted',30,angmin_def,
    c       See if this configuration is on the list
    c     print *,'ADJUST_XTRAJ ifirst,ilast,increment=',
    c    -  ifirst,ilast,increment
          character*200 trajnam12
          common /trajname/ trajnam12(2),ltrajnam12(2),ifirsttraj12(2),
          character*(*) label
          character*200 trajnam12
          common /trajname/ trajnam12(2),ltrajnam12(2),ifirsttraj12(2),
          character*8 lab12(2)
    c     print *,'WRITE_TRAJ_NAME incr=',incr,' iplot=',iplot,' iout=',iout
          character*30 talab(maxtors)
          character* 132 line(maxrec)
          character*1 ansrun
          character*4 atnami,atnamj
          character*8 resnam
          character*38 question
          call quiz(ansrun,iansrun,' ',' ',1,'torsion list input',18,0,
          call getresrange(nsegm,indexs,isegno,ixres,iresno,ifres,
    c       print *,'RESNO RANGE:',irefres1,irefres2
    c       print *,'ATOM RANGE:',ifres(irefres1),ilres(irefres2)
              call leftadjust4(line(index(ia))(ic1:ic2),atnami)
                  call leftadjust4(line(index(ja))(ic1:ic2),atnamj)
    c             Bond ia-ja found (on a heavy atom chain when nohrot=1)
                  call checktorbond(resnam,ixres(ia),ixres(ja),atnami,
    c             print *,'IA,ICAN,ICAC,IXRES=',ia,ican,icac,ixres(ia)
    c             Loop and normal torsion stepsizes
    c                 Bond ia-ja
                        call leftadjust4(line(index(iaa))(ic1:ic2),atnami)
                              call leftadjust4(line(index(jaa))(ic1:ic2),
            call getint('Number of torsions to track',27,0,1,49,ntang,0)
              call getintline(question,35,1,nslt,ixtor1234(1,it),4,0)
    c         print *,' ixtor:', (ixtor1234(k,it),k=1,4)
                call askyn('Do you want to use this angle',29,1,-1,iok,
            call blankout(talab(it),1,30)
          character*4 segid4(nmolslt)
    c     print *,'SETTORSLIM NSLT,NMOLSLT=',nslt,nmolslt
          call askyn('Do you want to set residue limits for torsions',46,
            call zeroiti(mask,0,nslt)
                call getint('Molecule number (0 to finish)',29,imol,1,
          character*1 ansrun
          character*4 extnam,lan(4),atnami,atnamj
          character*8 tname,resnam,rn
          character*200 inpfile,analfile
          character* 132 line(maxrec)
          character*2 iatnm2
          common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
          character*11 trajformatname
          common /trajectory/ nmmccheck,iftrajtyp(6),trajformatname(6)
          character*200 trajnam,trajnam2
          common /trajname/ trajnam,trajnam2,ltrajnam,ltrajnam2,ifirsttraj,
          character*4 namfcg
          character*4 tanames
          character*8 tnames
          common /tordat/ ntorn,tanames(4,28),tnames(28)
          common /connatdat/ ramax(99),ramax2(99),hlimfac,ianfg(99),
          common /graphics/ npixx,npixy,maxpixx,maxpixy,idwmain,idwplot,
          common /rotmat/ matrot0(4,4),matrot(4,4),nomat0
          common /columnlim/ incol(19),iidcol(19),iialtcol(19),iiinscol(19),
          character*5 crdext
          common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
          common /prokink/ icab(MAXHX),icaa(MAXHX),icb(MAXHX),ica(MAXHX),
    c     All arrays are of length maxframe
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
            call quiz(ansrun,nntyp,' ',' ',0,'Bond information source',23,0,
              call top_to_bond(nntyp,nneig,nhneig,ineig,iatnum,n,0,
              call nnlist(nslt,islvw,naslv,n,iatnum,ifchrg,c,nneig,nneiga,
            call quiz(ansrun,iansrun,' ','syntax for',10,
    c       if (extnam(2:3) .eq. 'mc') call askyn(
    c    -      'Do you want to include torsions moving hydrogens only',53,
    c    -      0,-1,nohrot,0,0)
            call openfile(40,0,'analysis',8,'new',analfile,namleno,
    c       MMC protein torsion input
    c         write (77,9781) ia,iatnum(ia),nneig(ia),nhneig(ia),mask(ia)
    c9781     format(i5,' iatno=',i2,' nn=',i2,' nnh=',i2,' mask=',i2)
                call leftadjust4(line(index(ia))(ic1:ic2),atnami)
    c             nhneig0=nhneig(ja)*nohrot
    c             if (ia .lt. ja .and. nneig(ia)-nhneig(ia) .gt. 1 .and.
    c    -            nneig(ja)-nhneig0 .gt. 1) then
                    call leftadjust4(line(index(ja))(ic1:ic2),atnamj)
    c               Bond ia-ja found (on a heavy atom chain when at most maxhrot
    c               hydrogens are moved)
                    call checktorbond(resnam,ixres(ia),ixres(ja),atnami,
    c             write (77,9877) ia,ixres(ia),resnam,ja,ixres(ja),atnamj,
    c    -           fixbond
    c9877         format(i5,1x,i5,1x,a,' - ',i5,1x,i5,1x,a,' FIXBOND=',i2)
    c               Loop and normal torsion stepsizes
                        call readint(line(index(ia)),iresncol1,iresncol2,
    c       Macromodel torsion input
    c             Bond ia-ja found
                          call leftadjust4(line(index(iaa))(ic1:ic2),lan(1))
                          call leftadjust4(line(index(ia))(ic1:ic2),lan(2))
                          call leftadjust4(line(index(ja))(ic1:ic2),lan(3))
                          call leftadjust4(line(index(jaa))(ic1:ic2),lan(4))
    c                       Protein, eliminate xxxxx, omega, and PRO phi,chi1,chi2
                              call leftadjustn(resnam,rn,8)
          character*1 ansrun
          character*25 read_format
          character*80 inpfiletmp,liner
    c     print *,'TOP_TO_BOND n=',n,' maxng=',maxng
            call quiz(ansrun,nntyp,' ',' ',0,'Bond information source',23,
    c       Get nnlist from a Charmm PSF file
                call openfile(inp_top,0,'PSF',3,'old',inpfiletmp,namlenn,
                call find_n_psf(inp_top,liner,nerr,n,natspsf,'!NATOM',6)
                    call askstop(1)
            call find_n_psf(inp_top,liner,nerr,0,natbond,'!NBOND',6)
            call checkdim(natbond,maxrec,'MAXREC',6,'number of bonds',15,0)
            call zeroiti(nneig,0,natspsf)
            call zeroiti(nhneig,0,natspsf)
    c           print *,'BOND:',i1,' - ',i2
    c         write (06,9782) i,iatnum(i),(ineig(j,i),j=1,nneig(i))
    c9782     format(i4,' iatno=',i2,' in=',10i4)
    c       Get nnlist from an Amber .top file
                call openfile(inp_top,0,'TOP',3,'old',inpfiletmp,namlenn,
            call find_ambertyp(inp_top,'%FLAG BONDS_INC_HYDROGEN',24,
            call zeroiti(nneig,0,n)
            call zeroiti(nhneig,0,n)
            call read_amber_bonds(1,inp_top,nneig,nhneig,ineig,iatnum,
              call find_ambertyp(inp_top,'%FLAG BONDS_WITHOUT_HYDROGEN',28,
            call read_amber_bonds(2,inp_top,nneig,nhneig,ineig,iatnum,
          character*25 read_format
          character*80 liner
    c     print *,'READ_AMBER_BONDS ih=',ih,' nr,ndig,maxng=',nr,ndig,maxng
          call blankout(liner,1,80)
            call blankout(liner,1,80)
              call lastchar(liner,lc,80)
                  call checkdim(nr,maxrec,'MAXREC',6,'number of bonds',
    c         write (6,8767) (i,itemp1(i),itemp2(i),i=nr0+1,nr)
    c8767     format(i5,' IT1,2=',2i6)
          character*4 segid4(nsegslt)
          character*8 resnames(maxrsd)
          character* 132 line(maxrec)
          character*(*) brslv,label
          character*4 atnam,bridgeats(200)
          character*1 ansrun
          character*4 bbhbats(6)
    c     print *,'GETHBANCHORDEF maxanchorlist=',maxanchorlist
    c     Establish solute atoms that are H-bond (chain) anchors
            call quiz(ansrun,iansrun,' ',label,llabel,
    c           Select backbone atoms
    c           Input atom name list
                call getnamelist(bridgeats,4,nbridgeats,'Anchor atom',11,
    c             print *,'ia,atnam=',ia,atnam
    c                 Anchor atom found
    c         Select by partial charge
              call definelist(ansrun,nslt,nanchor,ianchor,indexo,nsegslt,
            call askyn('Do you want to filter the list by charge',40,
            call getreal('Minimum (absolute) charge for an anchor atom',44,
            call zeroiti(indexa,0,nslt)
                call readreal(line(index(ia)),iiq1,iiq2,qa)
    c         Select from all solute atoms
    c         Reduce list
                call askyn('Do you want to repeat the selection',35,
          call getanchormod(ianchor2,iselfanc,nosameseg,iallsel,iw0)
    c           BB atom - skip
          call condenselist(ianchor,nanchor,0,iw0 )
    c     Remove carbons and aliphatic hydrogens
            call askyn('Do you want to repeat the selection',35,
          character*(*) label(6)
          character*4 segid4(mxrsd)
          character*132 line(mxrec)
    c     print *,'PRINTANCHORLIST ibondtyp=',ibondtyp,' iout=',iout
            call condenselist(ianchor,nanchor,0,iout)
            call condenselist(indexa,nanchorr,0,iout)
            call condenselist(indexov,nanchorn,0,iout)
          character*4 segid4(nsegslt)
          character* 132 line(maxrec)
          character*1 ansrun
          character*(*) bondname
    c     Establish solute atoms that are hydrophobic bond  anchors
    c     print *,'GETHPHANCHOR maxrec,maxanchorlist=',
    c    -         maxrec,maxanchorlist
          call zeroiti(indexa,0,nslt)
            call quiz(ansrun,iansrun,' ',' ',1,
            call definelist(ansrun,nslt,nanchorr,indexn,indexo,nsegslt,
            call getreal('Maximum (absolute) charge for an anchor atom',
    c     indexa will be nonzero for possible anchor atoms, negative if not anchor
    c       Hydrophobic 'bond'
              call askyn('Do you want to just use all carbons',35,1,+1,i0,0,
    c       Heavy-atom contact
    c       print *,'Anchor #',ja,'=',ia,' indexa(ia)=',indexa(ia),
    c    -    ' iatnum(ia)=',iatnum(ia)
                  call readreal(line(index(ia)),iiq1,iiq2,qa)
            call askyn('Do you want to repeat the selection',35,
          call getanchormod(ianchor2,iselfanc,nosameseg,iall,iw0)
          character*4 segid4(nsegslt)
          character* 132 line(maxrec)
          character*1 ansrun
          character*8 saltatname(100),atnam
    c     Establish solute atoms that are hydrophobic bond  anchors
    c     print *,'GETHPHANCHOR maxneig,maxrec,maxanchorlist=',
    c    -         maxneig,maxrec,maxanchorlist
          call zeroiti(indexa,0,nslt)
          call zeroit(chargesum,nslt)
          call quiz(ansrun,iansrun,' ',' ',1,
    c       Atom selection by charge
              call getreal('Minimum (absolute) charge for an anchor atom',
                  call readreal(line(index(ia)),iiq1,iiq2,chargesum(ia))
                  chargesum(ia)=charge(ia)
                    chargesum(ia)=chargesum(ia)+chargesum(ineig(in,ia))
    c       Atom selection by name
              call leftadjustn(atnam,atnam,8)
    c           Only oxygen is allowed to be negative
          call zeroiti(indexn,0,nslt)
            call quiz(ansrun,iansrun,' ',' ',1,
            call definelist(ansrun,nslt,nanchorr,indexn,indexo,nsegslt,
    c     indexa(ia) < 0: potential SB former
    c     indexn(ja): atom selected for anchor
            call askyn('Do you want to repeat the selection',35,
          call getanchormod(ianchor2,iselfanc,nosameseg,iall,iw0)
          character*4 segid4(nslt)
          character*1 ansrun
    c     print *,'GETMPXBDEF nslt,nsegslt=',nslt,nsegslt,' maxrsd=',maxrsd
          call zeroiti(indexa,0,nslt)
            call quiz(ansrun,iansrun,' ',' ',1,
            call definelist(ansrun,nslt,nanchorr,indexa,indexn,nsegslt,
            call quiz(ansrun,iansrun,' ',' ',1,
            call definelist(ansrun,nslt,nanchorn,indexov,indexn,nsegslt,
    c     Check for overlap
            call askstop(1)
          character*4 segid4(nsegslt)
          character*1 ansrun
          character*(*) label
    c     print *,'DEFINELIST nslt=',nslt,' NANCHORR=',nanchorr,
    c    -  ' MAXANCHORLIST=',maxanchorlist,' NSEGSLT=',nsegslt
              call indexit(indexn,1,nslt,0)
    c       Use all eligible solute atoms
            call indexit(indexn,1,nslt,0)
    c       Input anchor list
            call getlist(indexo,nanchoradd,1,nslt,1,maxanchorlist)
            call trnsfi(indexn(nanchorr+1),indexo,nanchoradd)
    c       Input segment number
            call getint('Segment number',14,1,1,nsegslt,iseganc,0)
            call indexit(indexn,nanchorr+1,nanchorr+nseg,
    c       Input anchor range
            call getrange(ifirst,1,ilast,nslt,increment,0,
            call indexit(indexn,nanchorr+1,nanchorr+nanchoradd,
    c       Input anchor residue list or range
              call getint('Segment number',14,1,1,nsegslt,iseganc,0)
    c       Input anchor residue list or range
    c         Get list
              call getlist(indexo,nanchorres,ifss,ilss,1,maxanchorlist)
    c         Get range
              call getrange(ifirstres,iresno(ifss),ilastres,iresno(ilss),
              call indexit(indexo,1,nanchorres,ifirstres-1)
              call findrange(iresno,ifss,ilss,iresanc,ifsr,ilsr,'residue',
    c           Add atoms ifsr - ilsr to the anchor list
    c           print *,'IR=',ir,' ADDING ',ifsr,' - ',ilsr
                call indexit(indexn,nanchorr+1,nanchorr+natadd,
    c           print *,'NANCHORR=',nanchorr
    c       print *,'Total number added to indexn=',nanchorr
            call askstop(0)
    c     print *,'DEFINELIST return'
          character*8 inout(2)
          character*9 onetwo(2)
            call askyn(
          call askyn('Do you want exclude intra segment/chain bonds',45,1,
          call getint('Minimum number of chemical bond to separate',43,
    c           Include (lev+1)th neighbors to the excluded list
    c               Check for duplicates
          character*4 atnami,atnamj
          character*8 resnami,resnam
    c     Check for small loops
            call leftadjustn(resnami,resnam,8)
    c     Get a list of atomic numbers that occured in this system
          call zeroiti(ixlist,0,99)
          call zeroiti(icatlist,0,15)
    c      write (6,7712) 'ialist',(ialist(i),i=1,nanos)
    c      write (6,7712) 'icatlist',(icatlist(i),i=1,nanos)
    c7712  format(1x,a,'=',15i3)
          character*2 iatnm2(99)
              call decidebondcut(ialist(i),ialist(j),rlim)
    c     Find the residue sequence number of residue number ir
    c     print *,'if,il=',ifirst,ilast
          character*(*) label
          common /logging/ logfile,ipredict
          character*80 line
    c     print *,'GETRESRANGE nslt,numres,nsegm,ipredict=',
    c    -  nslt,numres,nsegm,ipredict
              call reverseindex(indexs,iresno,ifres,1,numres,maxrec)
              call findrange(isegno,1,nslt,iseg1,ifss,ilss,'segment',7,
    c         print *,'ifss,ilss=',ifss,ilss
    c         print *,'iresno(ifss),iresno(ilss)=',iresno(ifss),iresno(ilss)
              call getint(line,lline,iresno(ifss),1,iresno(ilss),irn1,ihelp)
              call findrange(isegno,1,nslt,iseg2,ifss,ilss,'segment',7,
    c         print *,'ifss,ilss=',ifss,ilss
    c         print *,'iresno(ifss),iresno(ilss)=',iresno(ifss),iresno(ilss)
              call getint(line,lline,iresno(ilss),1,iresno(ilss),irn2,ihelp)
              call findsegres(isegno,iresno,ixres,1,nslt,iseg1,
              call findsegres(isegno,iresno,ixres,ia1+1,nslt,
            call askyn('Do you want to add an other range',33,1,-1,iadd,0,0)
          character*(*) lab
    c*****Find the representative atom of residue ir
          character* 132 line(maxrec)
          character*3 resnam3,repats,repat0
          character*4 atname,an
          character*8 resnam,rn
          common /represent/ repats(2,50),maxrepat
    c     Find first atom of residue ir
    c     write (77,*)'FINDAT START ir,ira1,ira2=',ir,ira1,ira2
    c     Find residue name and corresponding representative atom
          call leftadjustn(resnam,rn,8)
    c     print *,'resnam3=',resnam3,' resnam=',resnam
    c     Find representative atoms of residue ir
            call leftadjust4(atname,an)
    c     if (iarep .eq. 0 .and. infound .gt. 0) then
    c     do ia=ira1,ira2
    c       print *,'ia=',ia
    c       print *,line(index(ia))(1:78)
    c     end do
    c     stop
    c     end if
    c       Use all atoms
    c       Use the heavy atoms only
    c     Finds the PBC cell of the representative atom and shifts the whole molec
    c     back to the central cell
    c     print *,'PBCRESET ncell,ixyzexcld,ixyzincld=',
    c    -  ncell,ixyzexcld,ixyzincld
          call genimdist123dim(crep,cell,1,ncell,ixyzexcld,ixyzincld,
            call arrdiff(c(1,ia),cell(1,img),c(1,ia),3)
          call arrdiff(c2,c1,r,3)
          call genimdist(r,cell,1,ncell,img,d2)
    c*****Find the atom range in index containig ifind
          character*(*) label
    c     print *,'FINDRANGE n1,n,ifind=',n1,n,ifind
    c     print *,'FINDRANGE ifirst,ilast=',ifirst,ilast
    c*****Find the residue sequence number of segid=iseg, resid=ires
    c     print *,'FINDSEGRES ifst,ilst=',ifst,ilst
    c           print *,'FINDSEGRES iseg,ires,ia,irf=',iseg,ires,ia,irf
    c*****Allow for modification/addition of representative atomnames
          character*3 repats,resnam3,atnam3
          common /represent/ repats(2,50),maxrepat
          call askyn('Do you want to add/modify representative atom list',
            call getname(resnam3,len,
            call getname(atnam3,len,
          character* 132 line(maxrec)
          character*4 atnami,atnamj
    c     print *,'FINDPROTBACKBONE irespro,ia1=',irespro,ia1
            call leftadjust4(atnami,atnamj)
    c     print *,'ica,ic,in,icaonly=',ica,ic,in,icaonly
    c       Full backbone is required
    c       Only CA is required
    c         N and/or C is missing, CA is present - ask for CA only option
              call askyn('Do you want to look for just CAs',32,1,+1,icaonly,
          call zeroiti(ichiral,0,n)
              call listsame(ineig(1,ia),ian,4,n,nnpairs,nsamepair)
                  call comparetree(nneig,ineig,ian,ia,nnpairs(1,ipair),
    c     print *,'COMPARETREE ia,in1,in2=',ia,in1,in2
          call zeroiti(nian1,0,100)
          call zeroiti(nian2,0,100)
          call zeroiti(iused1,0,n)
          call zeroiti(iused2,0,n)
    c        write (6,7677) 'offsp1',(nnoffsp1(i),i=1,noffsp1)
    c        write (6,7677) 'offsp2',(nnoffsp2(i),i=1,noffsp2)
    c7677    format(1x,a,':',(30i3))
          character*2 iatnm2(99)
          character* 132 line(maxrec)
    c           Bond ia-ja found
    c                 iaa-ia-ja-jaa is a 1-4  chain
          common /nnwork/ ineig(MAXNEIG,MAXREC),nneig(MAXREC),
    c     print *,'RAMA init n,nres=',n,nres
          call zeroiti(nprossacc,0,36*nres)
          call zeroiti(issprossacc,0,5*nres)
          character* 132 line(mxrec)
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
          common /nnwork/ ineig(MAXNEIG,MAXREC),nneig(MAXREC),
          character*1 prosscode(5)
          character*4 atnam
    c     C - N - CA - C - N
    c     i1  i2  i3   i4  i5
    c         <--ires-->
    c     print *,'RAMA n,maxng,mxrec=',n,maxng,mxrec
          call trajlimtest(ixres(n),MAXFRAMES)
            call leftadjust4(atnam,atnam)
    c         Alpha carbon found
              call ca_to_bb(i3,iresno,nneig,ineig,index,line,ic1,
    c           Full residue backbone found
    c           write (iw0,*)'i1-5=',i1,i2,i3,i4,i5
    c           write (iw0,*)'ires i1-5=',iresno(i1),iresno(i2),
    c    -        iresno(i3),iresno(i4),iresno(i5)
    c               Save dial plot info
                        call trajlimtest(nframe,MAXFRAMES)
                call trajlimtest(nresfound,MAXFRAMES)
    c            write (iw0,2000) (c(k,i1),k=1,3),(c(k,i2),k=1,3),
    c     -       (c(k,i3),k=1,3),(c(k,i4),k=1,3),(c(k,i5),k=1,3)
    c2000         format(' c1=',3f8.3,' c2=',f8.3,' c3=',f8.3,' c4=',3f8.3,
    c     -         ' c5=',f8.3)
    c           Save and accumulate PROSS indices information
    c     Accumulate PROSS indices information
          call zeroiti(isspross,0,nresfound)
    c     Look for Alpha
    c     Look for Beta
    c     Look for turn
    c     Look for Pi
    c     Rest is Coil
    c1003  format(' WARNING: dial-plot residue #',i3,', residue',i6,' is ',
    c     -  'not a protein',/,10x,'No dials will be plotted for it')
          character* 132 line(mxrec)
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
          character* 132 line(mxrec)
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
    c     print *,'TORSIONDIALS n,mxrec=',n,mxrec
          call trajlimtest(nframe,MAXFRAMES)
          character* 132 line(mxrec)
          character*4 atnam
    c     Find the atom indices for C(n-1),N,CA,C,N(n+1)
    c     print *,'ic1,maxng,mxrec=',ic1,maxng,mxrec
    c     print *,'icaa,ires=',icaa,ires
              call leftadjust4(atnam,atnam)
    c     print *,' ina,ica=',ina,ica
    c       N-CA-C found
                call leftadjust4(atnam,atnam)
    c       print *,' icna   =',icna
                call leftadjust4(atnam,atnam)
    c       print *,' inca   =',inca
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
    c     print *,'RAMACHANDRANPLOT nres,ips,maxpres=',nres,ips,maxpres
    c         write (77,*) i,' phi,psi=',(res(k,i,maxpres),k=1,2)
            call rainbowscale(ips,ixm0,ixmm,25,nres,0.0,0.0,0.0,'N(res)',6)
            close (ips)
          character*(*) resnames(nres)
          common /nnwork/ ineig(MAXNEIG,MAXREC),nneig(MAXREC),
          character*1 prosscode(5)
          character*1 prossixy(36)
    c     print *,'RAMA hist nres=',nres
          character* 132 line(maxrec)
          character*2 iatnm2(99)
          character*1 sp
    c     print *,'BLS ic1,ic2,ir1,ir2=',ic1,ic2,ir1,ir2
    c     Bond length statistics first
    c         j<=i
    c     Bond angle statistics next
    c           k<=j
                cosa=dble(rjk)/sqrt(rj*rk)
          character* 132 line(maxrec),linep(25),blankline
          character*1 r1,bbsc(4)
          common /connatdat/ ramax(99),ramax2(99),hlimfac,ianfg(99),
          character*8 resnam
          character*8 atname,atname1,atnamea,atname1a
    c     print *,'HBLIST n,nslt,ic1,ic2=',n,nslt,ic1,ic2
    c     print *,'NHB=',nhbats
    c     print *,'hbf0,angm0=',hbf0,angm0
    c     print *,'hblimfac,angmin=',hblimfac,angmin
    c     print *,'isegno(1),(n),icallnn=',isegno(1),isegno(n),icallnn
    c     Determine if backbone,sidechain or other
            call changeprot(resnam,r1,2)
              call leftadjustn(atname,atnamea,8)
                call leftadjustn(atname1,atname1a,8)
    c          write (77,8871) i,atnamea,atname1a,indexn(i)
    c8871      format(i5,' atnamea=',a,' atname1a=',a,' indexn(i)=',i2)
          call checkhblist(n,ineig,nhbneig,maxng)
    c         Atom ia is always the donor H
    c         ibb,jbb: 1,2,3: Backbone/sidechain/solvent/?
    c         ihb0 is the heavy atom of the donor H
              call get_heavyat(ia,nneig,ineig,ixres,nframe,ihb0,maxng,
                  call angdistw(c(1,ia),c(1,ihb),c(1,ihb0),rHB,rb,rab,ang)
                  call readint(line(index(ia)),irn1,irn2,irna,2,1,irerr)
                  call readint(line(index(ihb)),irn1,irn2,irnb,2,1,irerr)
                call writeline(iwhbl,linep(ja),1,79,0)
    c         Cation - only H-bonds to water
    c           write (77,*) 'ia,ihb,nh12=',ia,ihb,nh12
                call angdistw(c(1,ihb),c(1,ia),c(1,nh12(1)),rHB,rb,rab,ang)
                  call angdistw(c(1,ihb),c(1,ia),c(1,nh12(2)),rHB,rb,rab,
                call readint(line(index(ia)),irn1,irn2,irna,2,1,irerr)
                call readint(line(index(ihb)),irn1,irn2,irnb,2,1,irerr)
                call writeline(iwhbl,linep(ja),1,79,0)
          character*132 line(n)
          character*22 bondname(5)
          character*(*) rlab
    c     print *,'HPHLIST n,maxng=',n,maxng
              call readint(line(index(ia)),irn1,irn2,irnia,2,1,irerr)
              call readint(line(index(ja)),irn1,irn2,irnja,2,1,irerr)
          character*132 line(mxrec)
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
          common /bondpairs/ ihbpair(2,MAXBONDS),ihb_pair_res(3,MAXBONDS)
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
    c     print *,'MPXBLIST n,iout=',n,iout
          call zeroiti(itemp,0,n)
          call zeroiti(itempres,0,n)
                    call readint(line(index(ia)),irnc1,irnc2,irnia,2,1,ire)
                    call readint(line(index(ja)),irnc1,irnc2,irnja,2,1,ire)
    c                 Save bond if not found yet
                          call askyn(
    c     Check for consistency of the neighbpor list
    c     Check for consistency of the H-bond list
    c     print *,'CHECKHBLIST n=',n
          character*(*) param,name
          character*8 resnames(maxrsd)
          character* 132 line(maxrec)
          common /nnwork/ sdm(MAX2D,MAX2D),nng(MAX2D),rmsd(MAX2D,MAX2D),
          call header_rrdist(iwrrdr,nrefrange,nnegrange,resdistlim,
          call header_rrdist(iwrrdc,nrefrange,nnegrange,resapplim,
          call header_rrdist(iwrrmp,nrefrange,nnegrange,0.0,
          call zeroiti(indexo,0,n)
          call zeroiti(indexn,0,n)
          call zeroiti(itemp1,0,n)
          call zeroiti(itemp2,0,n)
    c       First find the representative atom and atom range for residue irr
            call findat(irarep,ifres(irr),ilres(irr),line,index,
    c       write (iwrrdc,7711) irr,irarep,irra1,irra2,
    c    -     line(index(irarep))(ir1:ir1+3),
    c    -     line(index(irarep))(ic1:ic2)
    c7711  format('irr,irarep=',2i5,' irra1,irra2=',2i5,' r,a=',2a5)
                call findat(inarep,ifres(inr),ilres(inr),line,index,
    c           Calculate first based on representative atoms
                  call writeprox(iwrrdr,irarep,inarep,line,index,r2,
    c           Now obtain closest approach
                call findapproach(c,ifres(irr),ilres(irr),ifres(inr),
                  call writeprox(iwrrdc,irarepm,inarepm,line,index,
          call writeuniquelist(itemp1,ixresno,nres,resnames,ir2-ir1+1,
          call writeuniquelist(itemp1,ixresno,nres,resnames,ir2-ir1+1,
          call writeuniquelist(itemp2,ixresno,nres,resnames,ir2-ir1+1,
          call writeuniquelist(itemp2,ixresno,nres,resnames,ir2-ir1+1,
          call zeroiti(itemp1,0,n)
          call masktolist(itemp2,itemp1,nres,ncomplneg,0)
    c     itemp2 is the list of residues not on the neighbor and reference ist
                call findat(irarep,ifres(irr),ilres(irr),line,index,
                    call findat(inarep,ifres(irr1),ilres(irr1),line,
                      call findapproach(c,ifres(irr),ilres(irr),
                call findat(inarep,ifres(irr),ilres(irr),line,index,
                    call findat(irarep,ifres(irr1),ilres(irr1),line,
                      call findapproach(c,ifres(irr),ilres(irr),
    c     Calculate contact pairs
    c     Put the ref and neg atom lists into itemp3, itemp4. resp
    c     Find the proximal atoms for the ref and neg residues
          call zeroiti(indexo,0,n)
          call zeroiti(indexn,0,n)
    c           Mutually proximal atom pair found
                call writeprox(iwrrmp,ia,itemp1(ia),line,index,temp1(ia),
          call writeuniquelist(indexo,ixresno,nres,resnames,ir2-ir1+1,
          call writeuniquelist(indexn,ixresno,nres,resnames,ir2-ir1+1,
            call zeroiti(itemp,0,maxrsd)
            call condensemask(itemp,1,listrefres(nrefres),iw,maxrsd)
            call zeroiti(itemp,0,maxrsd)
            call condensemask(itemp,1,listnegres(nnegres),iw,maxrsd)
    c     print *,'MMDIST nmolslt,nmc,iw0,iw1=',nmolslt,nmc,iw0,iw1
            call zeroitd(c0,3)
                c0(k)=c0(k)+atw(ia)*c(k,ia)
                c0(k)=c0(k)+atw(ia)*c(k,ia)
              c2(k,im)=c0(k)/atwsum
          character* 132 line(maxrec)
    c      write (6,1000) index(irarep),index(inarep)
    c1000  format(' index(ir/n)=',2i4)
          character*(*) label
          character*8 resnames(mxrsd)
    c     print *,'WRITEUNIQUELIST ires1,ires2=',ires1,ires2,
    c    -  ' lab=',label(1:llabel)
          call zeroiti(itemp1,0,n)
          call zeroiti(itemp2,0,n)
    c     if (nmem .eq. 0) print *,'WRITEUNIQUELIST nmem=',nmem,
    c    -  ' ires1,ires2=',ires1,ires2
            call condensemask(index,ires1,ires2,iout,mxrsd)
            call condensemask(itemp2,ires1,ires2,iout,mxrsd)
          character*80 line
    c     print *,'CONDENSELIST len=',len
    c         Range found
              call printrange(line,i1,i2,ic,incr,iout)
          character*80 line
    c     print *,'CONDENSEMASK ires1,ires2=',ires1,ires2,' iout=',iout
    c      write (6,8888) (index(i),i=1,n)
    c8888  format(' index:',/,(10i5))
    c       Find first/next nonzero
    c         Find out if singleton or range
    c         print *,'i1=',i1,' i2=',i2,' ic=',ic
              call printrange(line,i1,i2,ic,0,iout)
          character*(*) line
    c       Singleton
            call writeint(line,ic,i1,lenw)
    c       Range
            call writeint(line,ic,i1,lenw)
            call writeint(line,ic,i2,lenw)
          character*8 resnames(mxrsd)
          character*(*) inpfile
          common /colorinfo/ ncolcode,maxcolcode
          common /nnwork/ sdm(MAX2D,MAX2D),nng(MAX2D),rmsd(MAX2D,MAX2D),
          character*4 yclab(1)
          character*14 distlab
          character*80 title2
          call getreal(title2,ltitle2,999999.0,rprtmax,1,0)
          call getreal(
          call getint('Number of colors to use',23,5,1,maxcolcode,ncolcode,
    c      do in=1,nresy
    c        write (iw,9855) in+inegres1-1,(rmsd(ir,in),ir=1,nresx)
    c9855    format(' in=',i4,200f6.2)
    c      end do
          call indexit(itemp,1,mxrsd,0)
          call plotmat(ips,kc,rmsd,dc,nresx,nresy,irefres1-1,inegres1-1,
    c     print *,'RCMIN,RCMAX=',rcmin,rcmax
          call colcodeminmax(ips,ixd,-60,0,ncolcode,maxcolcode,rcmin,rcmax)
          character*(*) angnames(ndials)
          character*(*) inpfile
          common /colorinfo/ ncolcode,maxcolcode
          common /nnwork/ ccc(MAX2D,MAX2D),itemp(MAX2D),fill(IFILL2)
          character*4 yclab(1)
          character*80 title2
          call getint('Number of colors to use',23,5,1,maxcolcode,
          call indexit(itemp,1,ndials,0)
          call plotmat(ips,kc,ccc,dc,ndials,ndials,0,0,0,0,1,0,0,iydel,00,
          call colcodeminmax(ips,ixd,-65,0,ncolcode,maxcolcode,rcmin,
          call pswrite(ips,150,25,'m',1)
          character*(*) angnames(ndials),label
          character*(*) inpfile
          common /nnwork/ ccc(MAX2D,MAX2D),irepav(MAX2D),irepmx(MAX2D),
          character*41 clstyp
          common /cluster_typ/ nclstyp,inumclst(9),ireadcutoff(9),
          character*1 ans
          character*80 line
              ccc(i,j)=1.0-abs(ccc(i,j))
            call getint('Number of clusters requested',28,999999,1,ndials,
            call getreal(line,22+llabel,999999.0,rdclust,1,0)
          call rmsdcluster(rdclust,1,ndials,index2d,iwt,ixclst,ifclst,
    c     Members of cluster ic: (index2d(i),i=ifclst(ic),ilclst(ic))
          call reportclust(ndials,0,1,nclust,ifclst,ilclst,index2d,value,
          character*8 resnames(mxrsd)
          common /nnwork/ rr1(MAX2D,MAX2D),nng(MAX2D),rr2(MAX2D,MAX2D),
          common /colorinfo/ ncolcode,maxcolcode
          character*4 yclab(1)
          character*200 plotfile,title2,inpfile1,inpfile2
          character*80 title1
    c     print *,'COMPARE_RRDIST maxcolcode,mxrsd=',maxcolcode,mxrsd
            call read_rrdist(rr1,41,'first',5,irefres1,irefres2,inegres1,
            call read_rrdist(rr2,41,'second',6,irefres11,irefres22,
          call openfile(42,0,'average distance matrix difference',34,'new',
          call askyn('Do you want absolute values in the difference plot',
          call getreal('Maximum distance to show',24,10000.0,dmaxshow,1,132)
          call zeroit(bfac,mxrsd)
          close (42)
          call openfile(ips,0,'average distance matrix difference',34,'new',
          call openps(ips,500.0,500.0,' ',1,' ',1,plotfile,0,
          call getint('Number of colors to use',23,5,1,maxcolcode,ncolcode,
          call indexit(itemp,1,mxrsd,0)
          call plotmat(ips,kc,rr1,dc,nresx,nresy,irefres1-1,inegres1-1,0,0,
          call colcodeminmax(ips,ixd,-60,0,ncolcode,maxcolcode,rdmin,rdmax)
          character*(*) lab
          character*80 inpmat,prompt,line
          call openfile(inpt,0,prompt,lprompt,'old',inpmat,linpmat,notfnd,
          close (inpt)
    cRes #     1(ILE ) - Res #     1(ILE ): <d>=  0.0000 A sd=  0.0000
          character*8 resnames(mxrsd)
          character(*) plotfile
          common /nnwork/ rnb1(MAX2D,MAX2D),nng(MAX2D),rnb2(MAX2D,MAX2D),
          common /colorinfo/ ncolcode,maxcolcode
          character*1 ans
          character*4 yclab(1)
          character*200 title1,title2,inpfile1,inpfile2
    c     print *,'COMPARE_BONDMAT maxcolcode=',maxcolcode,' MAX2D=',MAX2D
          call read_bondmat(rnb1,41,'first',5,nres,inpfile1,linpfile1,
    c     do iy=1,nres
    c       write (77,8811) iy,(rnb1(ix,iy),ix=1,nres)
    c     end do
          call read_bondmat(rnb2,41,'second',6,nres,inpfile2,linpfile2,
    c     do iy=1,nres
    c       write (78,8811) iy,(rnb2(ix,iy),ix=1,nres)
    c     end do
    c8811  format(' iy=',i5,/,(10f7.2))
          call openfile(42,0,'residue bond matrix difference',30,'new',
          call quiz(ans,idifftyp,' ',' ',0,'bond-difference plot type',25,
          call zeroit(bfac,nres)
    c       bfac(ir)=bfac(ir)/float(nres)
            call getreal('Minimum difference to show',26,rnbmin,rnbmin,1,0)
            call getreal('Maximum difference to show',26,rnbmax,rnbmax,1,0)
          close (42)
          call openfile(ips,0,'residue bond matrix difference',30,'new',
          call openps(ips,500.0,500.0,' ',1,' ',1,plotfile,0,
          call getint('Number of colors to use',23,5,1,maxcolcode,ncolcode,
          call indexit(itemp,1,mxrsd,0)
          call plotmat(ips,kc,rnb1,dc,nres,nres,0,0,0,0,1,0,0,iydel,00,
          call colcodeminmax(ips,ixd,-60,0,ncolcode,maxcolcode,rnbmin,
          character*(*) lab
          character*80 inpmat,prompt,line
          call openfile(inpt,0,prompt,lprompt,'old',inpmat,linpmat,notfnd,
          call indexit(itemp,1,nres,0)
    c     Read pointer if necessary
            call blankout(line,1,80)
    c       Read indices
            call lastchar(line,lc,80)
              call blankout(line,1,80)
              call lastchar(line,lc,80)
    c         Start of one column
    c           write (88,7211) iy,iyorig,ixf,ixl,
    c    -        (rnb(itemp(irx),iyorig),irx=ixf,ixl)
    c7211       format(' iy,iyorig=',2i5,' ixf,ixl=',2i5,' r=',10f7.2)
    c           write (88,7212) ixf,ixl,
    c    -        (rnb(itemp(irx),iyorig),irx=ixf,ixl)
    c7212       format(' ixf,ixl=',2i5,' r=',10f7.2)
          close (inpt)
          character*8 resnames(mxrsd)
          common /nnwork/ rmsf1(MAXRSD),rmsf2(MAXRSD),sd1(MAXRSD),
          common /colorinfo/ ncolcode,maxcolcode
          character*(*) plotfile
          common /t_test/ t_test_CI_1(5),t_test_CI_2(5),t_test_table(5,30)
          character*200 inpfile1,inpfile2
    c     print *,'COMPARE_RRDIST maxcolcode,mxrsd=',maxcolcode,mxrsd
          call read_rmsf(rmsf1,sd1,41,'first',5,irmin,irmax,
          call read_rmsf(rmsf2,sd2,42,'second',6,irmin2,irmax2,
          call openfile(42,0,'RMSF difference',15,'new',
          call zeroit(bfac,irmax)
    c     Both SD values are computed over a sample of 10 block averages
          character*(*) lab
          character*200 inpfile,prompt,line
          call openfile(inpt,0,prompt,lprompt,'old',inpfile,linpfile,notfnd,
          call zeroit(rmsf,mxres)
          call zeroit(sd,mxres)
          close (inpt)
    c    1 RMSF 16.96 SD=  1.99 decide=Uncorrelated
    c123456789012345678901234567890
    cRes #     1(ILE ) - Res #     1(ILE ): <d>=  0.0000 A sd=  0.0000
          character*1 marks(9)
          character* 132 line(maxrec)
          character*(*) inpfile,markfile
          common /nnwork/ rij1(MAXCONN,MAXCONN),rij2(MAXCONN,MAXCONN),
          character*8 fclab(10)
          character*200 title1,title2
    c     print *,'RRCONN numres1,numres2,nexpmax,iscalesum=',
    c    -  numres1,numres2,nexpmax,iscalesum
    c     print *,'RRCONN npspages,ipspage,maxrec=',
    c    -  npspages,ipspage,maxrec
    c       Distances based on representative atoms
    c         First find the representative atom and atom range for residue irr
              call findat(iarep(irr),ifres(irr),ilres(irr),line,index,
    c       Distances based on closest approach
                call findapproach(c,ifres(i1),ilres(i1),ifres(i2),
    c     Convert distances to adjacencies
    c     Calculate the adjaceny matrix powers
              colsum=0.0
                colsum=colsum+rij3(i1,i2)
            connmax=0.0
    c     Calculate correlation between marks and high/low col sums
          call rounddiv(numres,10,nx,nxdiv)
          call plotnps(xres,yres,MAXCONN,MAXCONN10,nexpplot,imf,iml,ifg,
          character*3 star3
          character*8 an1,an2,rn1,rn2,irn1,irn2
          character* 132 line(maxrec)
          character*80 pline
    c     For now, only maxddbin=20 works
    c     print *,'PAIRDISTPRINT nrescol,iresncol1,iresncol2=',
    c    -  nrescol,iresncol1,iresncol2
    c       dwavg=pairdistwsum(1,ip)/pairdistwsum(2,ip)
    c         if (maxddbin .le. 60) then
                  call blankout(pline,1,lastcol-1)
          character*(*) label
          call zeroiti(ihistogram,0,ng)
    c     imax=0
    c     do i=1,ng
    c       if (ihistogram(i) .gt. imax) imax=ihistogram(i)
    c     end do
    c     if (imax .gt. ng) then
    c       do i=1,ng
    c         ihistogram(i)=10.0*float(ihistogram(i))/float(imax)
    c       end do
    c     end if
          character* 132 line(maxrec)
          call bondcheck(iwchk,1,nslt,iatnum,nneig,ineig,maxng,c,maxdist,
    c       Mark segments to be disregarded
            call trnsfi(indexo,isegno,nslt)
          call contactcheck(iwchk,nslt,n,iatnum,c,line,irescol1,irescol2,
          character*9 decidebend
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
          character*60 message
          character*9 decide(6)
    c     print *,'HELIXAXIS maxhlx,nreshx,nrep=',maxhlx,nreshx,nrep
              calph(k,ir)=c(k,icaahx(ir))
              call dsmatvec(rot,calph(1,ir),calph(1,ir))
          call kahn(calph,nreshx,.true.,axisdir,axisini,axisend,rms,
              call calcperp(axisini,axisdir,calph(1,ir),camod(1,ir),
            call checkbend(calph,axisdir,camod,axfact,perpvec,helixcent,
            call calcturnperres(turnperres,nreshx,incrot,perpvec,axisdir,
            call trajlimtest(nframe,MAXFRAMES)
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
          character*2 ap_pa,in_ex
          common /signs/ itmem,normhx,isg2(MAXNHX2),memdir(MAXNHX),ap_pa(3),
    c     common /hxtrack/ ihxt1,ihxt2,torsav1(3,4),torsav2(3,4),x1(3),x2(3)
          character*80 line
    c     write (iw0,*) 'MULTIHELIX nhxres=',nhxres
          call blankout(line,1,80)
            cent1(1)=res(2,nframes,incr1+14)
            cent1(2)=res(1,nframes,incr1+15)
            cent1(3)=res(2,nframes,incr1+15)
              cent2(1)=res(2,nframes,incr2+14)
              cent2(2)=res(1,nframes,incr2+15)
              cent2(3)=res(2,nframes,incr2+15)
              call arrdiff(cent1,cent2,c12,3)
              c12_1=scprod(c12,ax1)
              c12_2=scprod(c12,ax2)
    c           Axes perpendicular
    c             Axes parallel
                  call project(cent2,cent1,ax1,ct1,a1)
              call paramx(cent1,ax1,a1,ct1)
              call paramx(cent2,ax2,a2,ct2)
              call arrdiff(ct1,ct2,ct12,3)
              c12norm=sqrt(scprod(ct12,ct12))
    c         write (iw0,8766) ihx,jhx,halflen1,halflen2,a1,a2
    c8766     format(' I,JHX=',2i2,' Hlen1,2=',2f6.2,' a1,2=',2f8.2)
    c           Make sure the closest approach points are inside the helix
                call paramx(cent2,ax2,-halflen2,s2)
                call paramx(cent2,ax2,halflen2,e2)
                call paramx(cent1,ax1,-halflen1,s1)
                call paramx(cent1,ax1,halflen1,e1)
                call project(s1,cent2,ax2,exs1,as1)
                call project(e1,cent2,ax2,exe1,ae1)
                call project(s2,cent1,ax1,exs2,as2)
                call project(e2,cent1,ax1,exe2,ae2)
    c             call compared2(s1,exs1,d2min,halflen2,-halflen1,a1,a2,
    c    -          as1)
    c             call compared2(e1,exe1,d2min,halflen2,halflen1,a1,a2,ae1)
    c             if (amin1(as1,ae1) .gt. halflen1) noproj=1
    c             call compared2(s2,exs2,d2min,halflen1,-halflen2,a2,a1,
    c    -          as2)
    c             call compared2(e2,exe2,d2min,halflen1,halflen2,a2,a1,ae2)
    c             if (amin1(as1,ae1) .gt. halflen2) noproj=1
    cn        call arrdiff(ct1,ct2,ct12,3)
    cn        ct12norm=sqrt(scprod(ct12,ct12))
    cn        do k=1,3
    cn          ct12(k)=ct12(k)/ct12norm
    cn        end do
                call trnsfr(torsats(1,4),c(1,icbahx(ihx)),3)
                call project(c(1,icbahx(ihx)),cent1,ax1,torsats(1,3),a)
                call arrsum(torsats(1,3),ax1,torsats(1,2),3)
    cn          call arrdiff(torsats(1,2),ct12,torsats(1,1),3)
                call project(cent2,cent1,ax1,s1,ac1)
                call arrdiff(cent2,s1,ct12,3)
                ct12norm=sqrt(scprod(ct12,ct12))
                  ct12(k)=ct12(k)/ct12norm
                call arrsum(torsats(1,2),ct12,torsats(1,1),3)
    c           write (iw0,3000)ihx,jhx,1,angrot1,((torsats(k,ii),k=1,3),
    c    -        ii=1,4)
    c           call trnsfr(torsav1,torsats,12)
                call trnsfr(torsats(1,4),c(1,icbahx(jhx)),3)
                call project(c(1,icbahx(jhx)),cent2,ax2,torsats(1,3),a)
                call arrsum(torsats(1,3),ax2,torsats(1,2),3)
    cn          call arrsum(torsats(1,2),ct12,torsats(1,1),3)
                call project(cent1,cent2,ax2,s2,ac1)
                call arrdiff(cent1,s2,ct12,3)
                ct12norm=sqrt(scprod(ct12,ct12))
                  ct12(k)=ct12(k)/ct12norm
                call arrsum(torsats(1,2),ct12,torsats(1,1),3)
    c           write (iw0,3000)ihx,jhx,2,angrot2,((torsats(k,ii),k=1,3),
    c    -        ii=1,4)
    c3000       format(' ihx,jhx=',2i3,' angrot',i1,'=',f8.1,' torsats=',
    c    -        4(2x,3f8.2))
    c           call trnsfr(torsav2,torsats,12)
    c         Helix distances
    c         End-end distances
              call dee(nframes,ixres-nhx*nhxres,ax1,ax2,cent1,cent2,
    c         Torsion angle over the closest aproach 'bond'
              call trnsfr(torsats(1,2),ct1,3)
              call trnsfr(torsats(1,3),ct2,3)
              call paramx(ct1,ax1,halflen1,torsats(1,1))
              call paramx(ct2,ax2,halflen2,torsats(1,4))
    c         write (77,6001) ihx,jhx
    c         write (77,6002) (i,' TA ',1,(torsats(k,i),k=1,3),i=1,4)
    c         write (77,6003) 'dhang',dhang
    c         Torsion angle over the center-center 'bond'
              call trnsfr(torsats(1,2),cent1,3)
              call trnsfr(torsats(1,3),cent2,3)
              call paramx(cent1,ax1,halflen1,torsats(1,1))
              call paramx(cent2,ax2,halflen2,torsats(1,4))
    c         write (77,6001) ihx,jhx
    c         write (77,6002) (i,' TA ',1,(torsats(k,i),k=1,3),i=1,4)
    c         write (77,6003) 'dhang_cc',dhang_cc
    c         if (ihx .eq. ihxt1 .and. jhx .eq. ihxt2) then
    c           if (nframes .eq. 1) write (78,6001) ihx,jhx
    c           write (78,6004) nframes
    c           call paramx(cent1,ax1,-halflen1,x1)
    c           call paramx(cent1,ax1,+halflen1,x2)
    c           write (78,6002) 1,' OB1',1,x1
    c           write (78,6002) 2,' OE1',1,x2
    c           write (78,6002) 3,' HC1',1,cent1
    c           write (78,6002) 4,' CA1',1,(c(k,iamin1),k=1,3)
    c           write (78,6002) 5,' CB1',1,(c(k,icbahx(ihx)),k=1,3)
    c           write (78,6002) 6,' N31',1,(torsav1(k,3),k=1,3)
    c           write (78,6002) 7,' N21',1,(torsav1(k,2),k=1,3)
    c           write (78,6002) 8,' N11',1,(torsav1(k,1),k=1,3)
    c           write (78,6002) 9,'CL1 ',1,ct1
    c           call paramx(cent2,ax2,-halflen2,x1)
    c           call paramx(cent2,ax2,+halflen2,x2)
    c           write (78,6002) 10,' OB2',2,x1
    c           write (78,6002) 11,' OE2',2,x2
    c           write (78,6002) 12,' HC2',2,cent2
    c           write (78,6002) 13,' CA2',2,(c(k,iamin2),k=1,3)
    c           write (78,6002) 14,' CB2',2,(c(k,icbahx(jhx)),k=1,3)
    c           write (78,6002) 15,' N32',2,(torsav2(k,3),k=1,3)
    c           write (78,6002) 16,' N22',2,(torsav2(k,2),k=1,3)
    c           write (78,6002) 17,' N12',2,(torsav2(k,1),k=1,3)
    c           write (78,6002) 18,'CL2 ',2,ct2
    c           write (78,6005) 1,2
    c           write (78,6005) 4,5
    c           write (78,6005) 5,6
    c           write (78,6005) 6,7
    c           write (78,6005) 7,8
    c           write (78,6005) 10,11
    c           write (78,6005) 13,14
    c           write (78,6005) 14,15
    c           write (78,6005) 15,16
    c           write (78,6005) 16,17
    c           write (78,6005) 9,18
    c           write (78,1000) 'ENDMDL'
    c         end if
    c1000 format(a)
    c6001 format('REMARK ihx,jhx=',2i3)
    c6002 format('ATOM  ',i5,1x,a4,1x,'TST',1x,'D',i4,1x,3x,3f8.3,
    c    -  '  1.0   0.0')
    c6003 format('REMARK ',a,'=',f10.4)
    c6004 format('MODEL',i6)
    c6005 format('CONECT',2i5)
    c     Find the nearest point p from x on c+a*ax
          character*2 ap_pa,in_ex
          common /signs/ itmem,normhx,isg2(MAXNHX2),memdir(MAXNHX),ap_pa(3),
          call paramx(cent1,ax1,-halflen1,s1)
          call paramx(cent1,ax1,+halflen1,h1)
    c       Establish 2nd axis direction wrt the first
            call paramx(cent2,ax2,-halflen2,s2)
            call paramx(cent2,ax2,+halflen2,h2)
          call paramx(cent2,ax2,-isg2(isg_i)*halflen2,s2)
          call paramx(cent2,ax2,+isg2(isg_i)*halflen2,h2)
          character*9 decidebend
    c     print *,'SOLUTEOVERLAY nsegslt,nslt,maxat=',nsegslt,nslt,maxat
          call zeroit(crmshift,3)
    c       No translation or rotation
            call trnsfr(cc2,c,3*nslt)
            call cofms(c,crmslt,nslt,atw)
            call arrdiff(crmslt,crmslt0,crmshift,3)
    c         Just translation
                call arrdiff(c(1,ia),crmshift,cc2(1,ia),3)
    c         Overlay of the whole solute
              call indexit(itemp,1,nslt,0)
              call bestoverlay(nslt,itemp,itemp,cres,c,atw,0.d0,
              call shiftmol(c,nslt,com2,cc2,-1.0)
              call rotate_c(cc2,nslt,rot,cc2,'HELIX',5)
              call shiftmol(cc2,nslt,com1,cc2,+1.0)
    c         Overlay separately each solute molecule
                call indexit(itemp,1,nats,0)
                call bestoverlay(nats,itemp,itemp,cres(1,ifat),c(1,ifat),
                call shiftmol(c(1,ifat),nats,com2,cc2(1,ifat),-1.0)
                call rotate_c(cc2(1,ifat),nats,rot,cc2(1,ifat),'HELIXs',6)
                call shiftmol(cc2(1,ifat),nats,com1,cc2(1,ifat),+1.0)
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
          character*2 iatnm2
          common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
          character*9 decidebend
    c     Compute the rotation angle of a helix form the average angle between
    c     the corresponding perpendiculars to the helix axis from the Calphas
    c     print *,'HELIXC nframe,nrep,iverbort=',nframe,nrep,iverbort
    c     cc2 will be the possibly translated/overlaid frame
          call helixaxis(cc2,nslt,0,calph,dirw,startw,endw,cent,perpvec,
          call trnsfrd(start,startw,3)
          call trnsfrd(end,endw,3)
          call trnsfrd(dir,dirw,3)
    c     Shift the helix so that the start is at the origin
            call dvdif(calph(1,ir),start,calph(1,ir))
          call dvdif(end,start,end)
          call zeroitd(start,3)
          call trnsfrd(startg,start,3)
          call trnsfrd(endg,end,3)
            call dvdif(cent,cent0,dev)
            call dvdif(startw,start0,devs)
            call dvdif(endw,end0,deve)
    c         For better result, rotate dir onto dir0. The rotation axis is the
    c         normal to the dir0-dir plane
              call dcross(dir0,dir,org)
                corig(k,1)=0.0
                ccurr(k,1)=0.0
                corig(k,2)=dir0(k)
                ccurr(k,2)=dir(k)
                corig(k,3)=org(k)
                ccurr(k,3)=org(k)
              call ormat(rot,ccurr,corig,3,iverbort,iw0)
                call dsmatvec(rot,calph(1,ir),calph(1,ir))
              call trnsfrd(dir,dir0,3)
              call calcperp(start,dir,calph(1,ir),org,perpvec(1,ir),it)
            changemin=100.0
            changemax=0.0
    c       write (79,*) 'NRES-',nres
              call angcomp(perpvec0(1,ir),dir0,perpvec(1,ir),angchange)
    c          write (79,9671) ir,angchange,(dir0(k),k=1,3),
    c     -      (perpvec(k,ir),k=1,3),(perpvec0(k,ir),k=1,3)
    c9671      format(' ir=',i3,' dang=',f10.4,' dir0=',3f10.5,/,
    c     -      ' perpvec=',3f10.5,' perpvec0=',3f10.5)
    cd77        write (77,6533) nframe,torsion*180.0/pi,changemin*180.0/pi,
    cd77     -    changemax*180.0/pi,(anglechange(ir)*180.0/pi,ir=1,nres)
    cd776533    format(' Nframe=',i6,' avg,min,max=',3f10.3,/,(10f8.2))
    c         Likely to have some sign flips
                call angcomp(perpvec0(1,ir),dir0,xx,angchange)
                  call trnsfrd(perpvec(1,ir),xx,3)
    c         Recalculate TPR
              call calcturnperres(turnperres,nres,incrot,perpvec,dir,
            call dcross(dir0,rn0,xx)
            call dcross(rn0,xx,yy)
    c       Keep the normal from oscillating 180 degrees
            call dsmatvec(rot,rn,xx)
            call printhelix(iw0,startw,endw,cent,rms,helixlen,dirw,angles,
    c     Calculate the angle between v1 and v2.
    c     Obtain the sign by requiring that d0 and d1 correspond to the z axis
          call dcross(v0,v,v01)
            call angcomp(perpvec(1,ir-1),axisdir,perpvec(1,ir),turnchange)
    c         Compare with turn angle in the reference conformation
    cd77          write (77,*) 'ir,turnchange,ref=',
    cd77     -      ir,turnchange*180.0/pi,anglechangeref(ir)*180.0/pi
    cd77        write (77,*) 'ir,turna,turnchange=',ir,
    cd77     -    turna*180.0/pi,turnchange*180.0/pi
    cd77      write (77,*) 'turnchangeav=',turnchangeav,turnchangeav*180.0/pi
          common /prokink/ icab(MAXHX),icaa(MAXHX),icb(MAXHX),ica(MAXHX),
          character*60 message
    c     ProKink variables
    c     print *,'PROKINKCALCLA icpr,inpr,nra,nrb=',icpr,inpr,nra,nrb
            calph(k,1)=c(k,icapr)
              calph(k,ir+1)=c(k,icaa(ir+1))
          call kahn(calph,nra+1,.true.,axisdira,axisinia,axisenda,
    c         calph(k,ir+1)=c(k,icab(ir+1))
              calph(k,ir)=c(k,icab(ir))
          call kahn(calph,nrb,.true.,axisdirb,axisinib,axisendb,rms,
            call fitpoints(prolinering,5,3,p0,ringnorm,axfact,iprintpk)
            call dvdif(pro_alphC,p0,xx)
    c     New code for proline kink calculation
    c     Bend angle: from the scalar product of the two axis vectors
    c     Both axis vectors point away from the proline
          cosa=-ddot(axisdira,axisdirb)
    c     Wobble  angle: from the scalar product of the two normals to the
    c     before helix axis (from the C-alpha of Proline and the end of the
    c     after helix)
          call calcperp(axisendb,axisdirb,pro_alphC,orig,perpvec,iprintpk)
          call dvsum(orig,axisdira,enda)
          call calcperp(axisinib,axisdirb,enda,origa,perpveca,iprintpk)
          call dvdif(origa,enda,xx)
          cc=dmag(xx)
          cosa=ddot(perpvec,perpveca)
    c     Establish sign
          call dcross(axisdirb,perpvec,zax)
          call calcperp(axisinib,axisdirb,bmin3C,orig3,perpvec3,iprintpk)
          call calcperp(axisinib,axisdirb,bmin4C,orig4,perpvec4,iprintpk)
          call dvnorm(perpvec)
          call dvnorm(perpvec3)
          call dvnorm(perpvec4)
          call dvnorm(perpvec34)
          cosa=ddot(perpvec,perpvec34)
    c     Establish sign
          cosa3=ddot(perpvec,perpvec3)
          cosa4=ddot(perpvec,perpvec4)
    c     For a line from start in the direction dir, calculate the normal to it
    c     from the point from. The normal meets the line at orig and its direction
    c     is perpvec
          call dvdif(from,start,dsx)
    c     call dcross(dsx,dir,xx)
    c     call dcross(xx,dir,perpvec)
          call dvnorm(perpvec)
              cc=ddot(dir,perpvec)
              call dvdif(start,orig,xx)
              cc=ddot(dir,xx)/dmag(xx)
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
          common /pcadat/ evecsprev(3,3),angprev(3)
          character*1 signlab(3)
    c     write (iout,*)'PRINCAX inputref,irefcall=',inputref,irefcall
          call trnsfr(co,c,3*n)
          call trnsfr(awo,aw,n)
          call extract(co,indxdel,3,n,nfinal)
          call extract(awo,indxdel,1,n,nfinal)
          call cofms(co,com,nfinal,aw)
          call zeroitd(tensinert,9)
    c     Find the eigenvectors a and eigenvalues mu of tensinert
          call dtred2(tensinert,3,3,diag,offdiag)
          call dtqli(diag,offdiag,3,3,tensinert,ierr)
          call indexit(index,1,3,0)
          call mrgsrt(6,index,evals,3,ifa,ila,itemp,temp,3)
    c       Test the eigenvectors after sort
    c       The columns of evecstemp should be also the eigenvectors
    c     The rows of the matrix evecs are the eigenvectors
    c      write (6,1000) diag,evals,tensinert,evecs
    c1000  format(' diag=',3f15.5,/,' evals=',3f15.5,/,
    c     -  ' tensinert=',/,3(3f10.5,/),' evecs=',/,3(3f10.5,/))
            call trnsfr(evals0,evals,3)
            call trnsfr(evecs0,evecs,9)
            call trnsfr(evecsprev,evecs0,9)
              call trnsfr(evals0,evals,3)
              call trnsfr(evecs0,evecs,9)
              call zeroit(ang,3)
          call indexit(index,1,3,0)
    c         Now check overlap with previous orientation
              call overlapcheck(evecsprev,evecs,overlap,index,nneg)
              call overlapcheck(evecs0,evecs,overlap,index,nneg)
            call overlapcheck(evecs0,evecs,overlap,index,nneg)
            call trajlimtest(nframe,MAXFRAMES)
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
          call zeroitd(rav,3)
          call zeroitd(rgten,9)
    c     Find the eigenvectors a and eigenvalues mu of rgten
          call dtred2(rgten,3,3,diag,offdiag)
          call dtqli(diag,offdiag,3,3,rgten,ierr)
            call trajlimtest(nframe,MAXFRAMES)
          call zeroit(overlap,9)
          call zeroiti(index,0,3)
          character*(*) name1,name2
    c       Calculate and print avg & sd
          corr=0.0
          character*(*) angname(maxdials),rname(maxrdat),title
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
          common /nnwork/ ccc(MAX2D,MAX2D),fill(IFILL2)
          character*26 corrtyp(2)
            cossum=0.d0
              cossum=cossum+res(1,i,incrres+id)
            cv(id)=1.0-sum/nframe
            cossum=cossum/sum
    c         Calculate circular correlation coefficient between angles
                ccc(id,id)=1.0
    c               Fisher & Lee (slow O(nframe^2))
    c               Jammalamadaka and SenGupta (2001) (fast, O(nframe))
                  ccorr=(sinsum_ab)/dsqrt(sinsum2_a*sinsum2_b)
                  ccc(id,jd)=ccorr
                  ccc(jd,id)=ccorr
    c         Calculate standard correlation coefficient between regular data
                  call correl(iout,res(1,1,incrres+id),ic,
    c         Calculate linear-circular correlation coefficient between
    c         angles and regular data (Mardia)
                  corrsum=0.d0
                  cossum=0.d0
                  cossum2=0.d0
                  crc=0.d0
                  crs=0.d0
                  csc=0.d0
                    c=res(1,i,incrres+id)
                    crc=crc+x*c
                    crs=crs+x*s
                    csc=csc+s*c
                    cossum=cossum+c
                    cossum2=cossum2+c**2
                    corrsum=corrsum+ang_ai*x
                  corrsum=corrsum*radtodeg
                  cav=cossum/nframe
                  csd=(cossum2/nframe-cav**2)
                    csd=dsqrt(csd)
                    ccorr=(r12**2+r13**2-2.0*r12*r13*r23)/(1-r23**2)
                    ccorr=0.0
    c     write (77,*) 'ndials=',ndials
    c     do i=1,ndials
    c       write (77,1004) (ccc(i,j),j=1,ndials)
    c     end do
    c1004 format(5e13.6)
    c     Perform the Monte Carlo estimate of ligand and interfacial volumes
          character*80 line
          call cellpart(c,ian,itemp1,1,na,padding,spacing,corner,ecell,edge,
    c     Get a list of filled cells
    c     Mark neighbors of filled cells
          call zeroiti(itemp1,0,ntotcell)
            call randpx(3,rand)
    c         Shift random copy to cell icell and examine just its neighbor cells
    c         Deconvolute icell
              call zeroiti(iamn,0,nsegslt)
    c                 Cell is not empty
    c           Check if potential interface neighbor
    c           Inside solute's vdW region
    c           Within solvent-excluded shell
                  call blankout(line,1,80)
    c             write (77,2000) line(1:80)
    c           Within first solvation shell
    c           Possible interface
    c                   Calculate the scalar product
    c                   Calculate the scalar threshold value
                    call blankout(line,1,80)
          character*8 resnames(mxrsd)
          character*1 minmaxlab1(mxrsd)
          character*1000 line
          character*30 title(MAXCOL),label(4),syslabel
          character*30 source(3)
          character*1 ans,types1(4)
          character*26 err_typ(2)
          character*80 atitle
          character*200 inpstatfile,outtabfile
          common /nnwork/ data(MAXCOL,MAXREC),err(MAXCOL,MAXREC),
    c     print *,'SUMMARIZE_AMBER_CSV inpt=',inpt
          call blankout(inpstatfile,1,80)
            call openfile(inpt,0,'statistics',10,'old',inpstatfile,
          call openfile(iout,0,' ',1,'new',outtabfile,
          call blankout(line,1,len)
            call blankout(line,1,len)
          call blankout(atitle,1,80)
          call lastchar(atitle,latitle,80)
          call getskipcomma(inpt,line,len,llab,ifail)
          call findname(syslabel,label,1,4,nc,lsyslabel)
              call getskipcomma(inpt,line,len,llab,ifail)
          call quiz(ans,ns,'t',' ',0,'Data source to tabulate',23,0,5,6,0)
            call getskipcomma(inpt,line,len,llab,ifail)
    c     lenergy_source=llab
    c     energy_source(1:lenergy_source)=line(1:llab)
          call blankout(line,1,len)
          call findnextchar(',',line,ic,len)
          call findnextchar(',',line,ic,len)
            call findnextchar(',',line,ic,len)
    c       print *,'IC=',ic
              call blankout(title(ncol),1,30)
          call askyn('Do you want to tabulate the errors too',38,1,-1,
          call askyn('Do you want to mark with m and M the extreme values',
          call read_amb_data_csv(data,err,ncol,ixr1,ixr2,resnames,nresdat,
              call pickname('Index of choice (0 to finish the list)',38,
                  call askyn('OK',2,1,-1,iok,0,0)
                  call askyn('OK',2,1,+1,iok,0,0)
          call zeroiti(iresix,0,maxresno)
          call zeroit(colsumm,ncolp)
    c         Establish min and max
                colsumm(ip)=colsumm(ip)+data(ixprint(ip),ir)
    c       Pairs
              call getrange(ires1,ixr1(1),ires2,maxresno,incr,0,
              call zeroit(colsum,ncolp)
    c         write (6,8792) (ixprint(ip),ip=1,ncolp)
    c8792     format(/,' IXPRINT=',10i4)
    c             Requested range found
    c               Obtain the min/max labels for this residue
                  call zeroit(colsum,ncolp)
                        colsum(ip)=colsum(ip)+data(ixprint(ip),irr)
                      colsumm(ip)=colsumm(ip)+colsum(ip)
                        colsum(ip)=colsum(ip)+data(ixprint(ip),irr)
                      colsumm(ip)=colsumm(ip)+colsum(ip)
    c         Single property - print residue-residue matrix
                call getrange(iresrow1,ixr1(1),iresrow2,maxresno,incr,0,
                call restoix(iresix,iresrow1,iresrow2,ix1,ix2,ifail12,
                call getrange(irescol1,ixr1(1),irescol2,maxresno,incr,0,
                call restoix(iresix,iresrow1,iresrow2,ix1,ix2,ifail,
                call askyn('Do you want to change the ranges',32,1,-1,ichng,
    c         write (6,*) 'irescol1,irescol2,iresrow1,iresrow2=',
    c    -             irescol1,irescol2,iresrow1,iresrow2
    c         write (6,*) 'irescol1p,irescol2p,iresrow1p,iresrow2p=',
    c    -             irescol1p,irescol2p,iresrow1p,iresrow2p
    c         write (6,*) 'nres,ntabcol,ntabrow=',nres,ntabcol,ntabrow
                call askyn('Do you want to change the ranges',32,1,1,ichng,
    c          write (6,7711) 'ixr1:',(ixr1(ii),ii=1,nres)
    c          write (6,7711) 'ixr2:',(ixr2(ii),ii=1,nres)
    c          write (6,7711) 'ixrow:',(ixrow(ii),ii=1,ntabrow)
    c          write (6,7711) 'ixcol:',(ixcol(ii),ii=1,ntabcol)
    c7711      format(1x,a,/(20i4))
              call zeroit(colsumm,ntabcol)
              call zeroiti(ixcolmin,0,ntabcol)
              call zeroiti(ixcolmax,0,ntabcol)
              call zeroit(tabrowsum,ntabrow)
                colmin(irc)=999999.0
                colmax(irc)=-colmin(irc)
                  colsum(irc)=colsum(irc)+r
                    colmin(irc)=r
                    colmax(irc)=r
    c             write (77,8972) irr,if0,ixrow(irr),irc,ixcol(irc),r
    c8972         format(' irr=',i4,' if0=',i5,' ixrow(irr)=',i5,' irc=',i4,
    c    -          ' ixcol(irc)=',i4,' data=',f10.6)
              call askyn('Do you want to tabulate an other energy term',44,
          call askyn(
    c       print *,'MAXRESSNO=',maxresno,' NRES=',nres,' IPAIRS=',ipairs
    c       print *,'IXR1:'
    c       write (6,8768) (ixr1(ir),ir=1,nres)
    c       print *,'IXR2'
    c       write (6,8768) (ixr2(ir),ir=1,nres)
    c       print *,'IRESIX:'
    c       write (6,8768) (iresix(ir),ir=1,maxresno)
    c8768   format(i3,19i4)
              call pickname('First property to correlate (zero to exit)',42,
                call pickname('Second property to correlate',28,title,
                  corr=(s12-s1*s2/nres)/
                    call getint('Residue number of the selected residue',38,
    c             Find the row numbers corresponding to the resid numbers irr1,irr2
                    col1(i-ixfres(ir1)+1)=data(ic1,i)
                    col2(i-ixfres(ir1)+1)=data(ic2,i)
                    corr=(s12-s1*s2/nrows)/
            call askyn(
          call blankout(atitle,1,80)
          call lastchar(atitle,latitle,80)
          call getskipcomma(inpt,line,len,llab,ifail)
          call findname(syslabel,label,1,4,nc,lsyslabel)
    c2003  format(' A PDB file can be also read to provide residue names',/,
    c     -  ' Residue numbers should correspond to the residue numbers ',
    c     -  'read by mmpbsa.pl',/,
    c     -  ' Hitting enter will skip reading/using residue names')
          character*8 resnames(mxrsd)
          character*3 rn1,rn2
          character*1000 line
    c     Read data from Amber energy analysis tables
    c     print *,'READ_AMB_DATA_CSV ipairs,ncol=',ipairs,ncol
            call blankout(line,1,len)
            call lastchar(line,lc,len)
    c         write (6,*) 'NRES=',nres,' LINE=',line(1:20)
              call findnextchar(',',line,ic,len)
                call findnextchar(',',line,ic,len)
                call blankout(resnames(ix1(nres)),1,8)
                call findnextchar(',',line,ic,len)
                call blankout(resnames(ix1(nres)),1,8)
                call blankout(resnames(ix2(nres)),1,8)
    c           write (77,8977) nres,rn1,rn2,ix1(nres),ix2(nres),
    c    -       resnames(ix1(nres)),resnames(ix2(nres))
    c8977        format(' nres=',i4,' rn1,2=',a,'|',a,'| ix1,2=',2i5,
    c    -         ' rn1,2=',a,'|',a,'|')
                call getnextcsv(line,i,data(nc,nres),2,ic,len)
                call getnextcsv(line,i,sd,2,ic,len)
                call getnextcsv(line,i,err_mean,2,ic,len)
    c         write (77,8877) ncol,nc,nres,(data(i,nres),i=1,ncol)
    c8877     format(' ncol,nc,nres=',i2,i3,i5,' data=',10f10.3)
          character*1000 line
          character*10 label,labelread
          call blankout(line,1,80)
          call lastchar(line,ilc,80)
    c     write (77,*) 'CHECKLABEL',line(1:ilc)
    c     print *,'RESTOIX ir1,ir2=',ir1,ir2,' ix1,ix2=',ix1,ix2
          character*1000 line
            call blankout(line,1,len)
            call lastchar(line,lc,len)
    c     print *,'GETSKIPCOMMA lc=',lc,' line:',line(1:lc)
          character*1 minmaxlab1(mxrsd)
    c     print *,'LABMINMAX ir,ncol=',ir,ncol
          character*(*) line
          call nextchar(line,ic,len)
          call nextblank(line,ic,len)
          character*(*) line
          call nextchar(line,ic,len)
          call nextblank(line,ic,len)
    c     print *,'GETNEXTR ic1,ic2=',ic1,ic2
          character*(*) line
          call findnextchar(',',line,ic,len)
    c       print *,'GETNEXTR ic1,ic2=',ic1,ic2
    c       Last value, without comma
            call nextblank(line,ic,len)
          character*(*) q,lab,list(nlist)
    c     print *,'FINDINDEX llab=',llab
    c     Partition the solute atoms into cells for fast distance calculation
    c     print *,'CELLPART NA1,NA,SPACING,PADDING=',na1,na,spacing,padding
          call extension(c,ian,0,na1,na,xyzmin,xyzmax,cent,0,0,v)
            corner(k)=xyzmin(k)-padding
    c     Establish the cell of each solute atom
          call mrgsrti(6,index,icellno,nha,ifirst,ilast,itemp1,itemp2,
          call zeroiti(ifirst,0,ntotcell)
    c     Establish the limits of cells
    c     For empty cells, ifirst(icell)=0!
    c     print *,'return CELLPART NA=',na
          character*(*) labdial(ndials)
          character*80 title
          character*(*) remark
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
          character*200 trajnam,trajnam2
          common /trajname/ trajnam,trajnam2,ltrajnam,ltrajnam2,ifirsttraj,
    c     print *,'DIALPS idincr=',idincr,' ndials=',ndials
    c     xm=500.0
            call openps(iout,xm,ym,title(1:ltitle),ltitle,remark,lremark,
            call plothead(iout,xm,ym,title,ltitle,remark,lremark,trajnam,
            call partwindow(xm,ym,ymfac,nddo,edge,x0,y0,nddone,ndprow,nrows,
              call drawdial(iout,edge,nframesaved,nfravgd,id,labdial(id),
    c     The window/page is xm by ym, (0,0) is the bottom left corner
    c     The subroutine finds the upper left corners (x0(i),y0(i)) of nd cubes
    c     that best fill the window/page
    c     print *,'PARTWIN nd,ndprow=',nd,ndprow
    c       x00=(ndprow-ndo)*(xm-ndo*edge)/2.0
    c     write (6,1000) xm,ym,edge,(i,x0(i),y0(i),i=1,nd)
    c1000  format(' xm=',f10.4,' ym=',f10.4,' edge=',f10.5,/,
    c     -  i4,' x0=',f10.4,' y0=',f10.4)
          character*(*) lab
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
    c     Draw a 'dial' into the square of edge edge, with upper left corner at
    c     (x0,y0), title lab(1:llab)
    c     print *,'DRAWDIAL n=',n,' incrid,ll=',incrid,llab,' idial=',idial
    c     print *,'DRAWDIAL ioutpr,mappdf=',ioutpr,mappdf
          cx=x0(idial)+edge*0.50
    c     if (ndprow .gt. 4) ldown=8
          cy=y0(idial)-edge*0.50+lyshift
    c     Print 0,90,180,270 marks
            call psshow(iout,'0',1)
          call psshow(iout,'p/2',3)
            call psshow(iout,'p',1)
          call psshow(iout,'3p/2',4)
          call psshow(iout,lab,llab)
    c     Draw the circle
    c     Draw ticks
    c     Draw axes
    c     Draw inner disk and initial postion line
    c     write (iout,*) '0.7 0.7 0.7 setrgbcolor'
          call rgbcolor(iout,-4)
          cosav=0.d0
    c     write (77,*) 'DRAWDIAL: ',lab(1:llab)
          call rgbcolor(iout,-6)
    c         First average over incrloop steps
              cosavinc=0.0
                cosavinc=cosavinc+res(1,j,incrid+idial)
              cosavinc=cosavinc/incrloop
              cosavinc=cosavinc/sqsum
    c           if (angdev .gt. pi) angdev=angdev-2.0*pi
    c           if (angdev .lt. -pi) angdev=angdev+2.0*pi
    c           write (77,8788) ang*rdtodg,angprev*rdtodg,angdev*rdtodg,
    c    -        angdevrev*rdtodg,nd
    c8788       format(' ang,prev=',2f7.1,' angdev,rev=',2f7.1,' nd=',i4)
    c             if (lab(1:1) .eq. 'P') print *,'i,ri,angdev,nd=',
    c    -            i,ri,angdev,nd
    c         Just draw dots
            cosav=cosav+res(1,i,incrid+idial)
          cv=1.d0-dsqrt(cosav**2+sinav**2)/dfloat(n)
    c       Draw average line
            call rgbcolor(iout,-1)
            cosav=cosav/n
            cosav=cosav/sqsum
    c       Draw tick outside at the final position
            call rgbcolor(iout,-6)
    c     Draw CV bar
          call rgbcolor(iout,-9)
          call rgbcolor(iout,9)
            call zeroiti(npdf,0,360)
            call zeroit(pdf,360)
    c         ix=(ang-angmind)/grd+1
    c         if (ix .gt. 360) ix=360
    c       ixmax=(angmaxd-angmind)/grd+1
    c       if (ixmax .gt. 360) ixmax=360
    c       write (88,*) 'R1=',r1,' ANGMIND,MAXD=',angmind,angmaxd
    c       write (88,*) 'ANGMIN,MAX=',angmin,angmax
    c           write (88,*) ix,' ANG=',ang
          character*(*) string
          character*200 clean
              clean(ic+inc:ic+inc)='\\'
              clean(ic+inc+1:ic+inc+1)='('
              clean(ic+inc:ic+inc)='\\'
              clean(ic+inc+1:ic+inc+1)=')'
              clean(ic+inc:ic+inc)=string(ic:ic)
          character*(*) title,title1,filename,filename2
    c     print *,'OPENPS'
          call psheader(iout,title,ltitle,-10,-10,ixm,iym,npspages,ipspage)
          call plothead(iout,xm,ym,title,ltitle,title1,ltitle1,
          character*(*) title
          character*12 today
          common /today_date/ ltoday,today
    c     print *,'PSHEADER npspages,ipspage=',npspages,ipspage
    c       if (npspages .lt. 10) then
    c         write (iout,1021) npspages
    c       else if (npspages .lt. 10) then
    c         write (iout,1022) npspages
    c       else
    c         write (iout,1023) npspages
    c       end if
    c1021  format('%%Pages:',i1)
    c1022  format('%%Pages:',i2)
    c1023  format('%%Pages:',i3)
          character*(*) title,title1,filename,filename2
            call psshow(iout,title1,ltitle1)
              call psshow(iout,'Simulaid-generated plot',23)
            call lastchar(title,lct,ltitle)
            call psshow(iout,title,lct)
          character*8 resnames
          character*8 brslv
          character*80 line
    c     write (iout,*) 'ianchor2,iselfanc,n,ixres(1),ixres(n)=',
    c    -  ianchor2,iselfanc,n,ixres(1),ixres(n)
    c     print *,'HBBRIDGE nanchor=',nanchor,' maxrec=',maxrec
    c       write (iout,2311) iat,indexa(iat),resnames(ixres(iat)),
    c    -    (ineig(maxng+1-ia,iat),ia=1,nhbneig(iat))
    c2311    format(' HBBR ianch=',i6,' indexa=',i4,' resn=',a,' ihb=',5i6)
    c       Grow tree at depth maxbrlen (4)
    c           ixp is the parent residue (solvent molecule) number
    c             Drop if it is a loop back or already in the current list
    c             list is only sorted at each level
    c               print *,'it,itt,iflist(itt),illist(itt)=',
    c    -             it,itt,iflist(itt),illist(itt)
                    call findixsort(list,iflist(itt),illist(itt),ian,
    c             do j=1,il-1
    c               if (nnblist(j) .eq. ian) go to 100
    c             end do
                  call findixsort(nnblist,1,nn,ian,ixian,itrynnb)
    c             If anchor atom then finish bridge and gather statistics
    c               See if end is also an anchor
                      call findixsort(ibridgetype(1,it,iia),1,
    c                 do ial=1,nbridgetype(it,iia)
    c                   if (ibridgetype(ial,it,iia) .eq. ian) then
    c                     iend=ial
    c                     go to 300
    c                   end if
    c                 end do
    c                   Not found, add to list
    c                   Shift ibridgetype,lpath to maintain them sorted
    c                   ibridgetype(nbridgetype(it,iia),it,iia)=ian
    c                   iend=nbridgetype(it,iia)
    c                   write (iout,1002) ibridgetype(iend,it,iia)
    c                     write (iout,*) 'TBCK it,iit,ixp,iparent(ixp)=',
    c    -                 it,iit,ixptb,iparent(ixptb)
                              call blankout(line,1,ic)
    c                  write (77,8877) il,ia,iend,it,iia,lpath(iend,it,iia)
    c8877              format(' il,ia,iend,it,iia=',5i4,' lpath=',i3)
    c               New atom, add to list of next level all atoms of this solvent
    c               write (iout,*) 'ADD ian,isolv,iparent(isolv)=',
    c    -                              ian,isolv,iparent(isolv)
    c               Add all atoms of this residue to nnblist
    c               write (iout,*) 'Adding itrynnb,ifr,ilr,nn=',
    c    -                          itrynnb,ifr,ilr,nn
    c                 do isv=ifr,ilr
    c                   nnblist(nn+isv-ifr+1)=isv
    c                 end do
    c               nn=nn+ilr-ifr+1
    c          if (nn .gt. 0) write (iout,6633) (nnblist(kk),kk=1,nn)
    c6633      format(' nnblist=',10i6)
              call trnsfi(list(ll+1),nnblist,nn)
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
          character*80 title
          character*200 pdbout
    c     Print a PDB file of atoms where each cluster is a separate residue
    c     print *,'PRINTCLUSTERPDB ncl=',ncl,' na=',na
          call openfile(50,0,'clustered PDB',13,'new',pdbout,nl,notfnd,
          close (50)
          character*(*) outfile
          character*1 xyzl
          character*200 centfile
          common /axislab/ xyzl(3)
    c     print *,'PRINTCLUSTEREXT ncl=',ncl,' na=',na,
    c    -  ' iclstyp,iout=',iclstyp,iout
    c         Convert xyz to polar coordinates
              call polar(dxyz(1),dxyz(2),dxyz(3),r,phi,theta)
    c         Calculate cluster COM
              call zeroit(cent(1,icl),3)
                  cent(k,icl)=cent(k,icl)+c(k,ixclst(ia))
                cent(k,icl)=cent(k,icl)/float(ilclst(icl)-ifclst(icl)+1)
          call askyn('Do you want a cluster-center file',33,
            centfile(1:namleno)=outfile(1:namleno)
            centfile(namleno+1:namleno+4)='.cnt'
            call openfile(51,0,'cluster center',14,'new',centfile,namleno+4,
          character*(*) bondname
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
          common /bondpairs/ ihbpair(2,MAXBONDS),ihb_pair_res(3,MAXBONDS)
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
    c     dimension ixx(3000)
    c     ianchor: list of anchor atoms
    c     indexa(ia): one if ia is an anchor atom;
    c                 -1 if it is non-anchor atom but allowed by the charge filter
    c     iqfsel2: If > 0, apply the charge filter to both atoms in the H bond
    c     ianchor2: If > 0, both atoms forming the H bond have to be anchors
    c     iselfanc: If = 0, exclude anchor-anchor hydrogen bonds
    c     write (iout,*) 'SELECTBOND nframe,nbresfound,ianchor2,iselfanc=',
    c    -  nframe,nbresfound,ianchor2,iselfanc
          call trajlimtest(nframe,MAXFRAMES)
          call zeroiti(itemp,0,mxbonds)
          call zeroiti(itempres,0,mxbonds)
    c         See if end is also an anchor
    c           write (iout,7943) ian,iat,indexa(ian),ib1,ib2
    c7943       format(' ian,iat=',2i5,' indexa(ian)=',i4,' ib1,2=',2i5)
                    call askyn('Do you want to continue without tracking',
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
    c     print *,'TRACKSTAT nbonds=',nbonds,' nframe=',nframe
          call zeroiti(nhbdist,0,nbonds)
          call zeroiti(nhbpers,0,nbonds)
          call zeroiti(maxlenon,0,nbonds)
          call zeroiti(maxlenoff,0,nbonds)
          call zeroiti(isegstart,0,nbonds)
          call zeroiti(itf,0,nbonds)
          call zeroiti(itl,0,nbonds)
            call readbitc(ires(1,ifr),it,nbonds,30,MAXITEMS)
            call trnsfi(itprev,it,nbonds)
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
    c     Read bits
    c     print *,'GETBONDTRACK itr,nbits,nfrm,nfram=',itr,nbits,nfrm,nframe
          call zeroiti(itrack,0,nframe)
          common /aucw/ auc_all(MAXFRAMES21,MAXCOPY)
          character*6 frunit(4)
    c     Calculate the autocorrelation of bond itr
    c     write (6,*) 'AUTOCORR itr_ix,itr,ifirstframe,lastframe=',
    c    -  itr_ix,itr,ifirstframe,lastframe
          call zeroit(auc,nframe2)
          call zeroiti(nauc,0,nframe2)
    c      write (0006,9867) itr,nframe,lastframe,loffmin
    c9867  format(' ITR=',I6,' NFRAME=',I6,' LASTFR=',I6,' LOFFMIN=',I6)
    c     Accumulate auc
              call zeroiti(itrack,nframe,lastframeinp)
    c         Pad with repeating the track
    c          write (iout,*) 'NADD=',nadd
    c            write (iout,9682) incr,incr0,nadd,ncop
    c9682        format(' INCR=',I6,' INCR0=',I6,' NADD=',I6,' NCOP=',I6)
    c           print *,'NCOP,NSEG_SCR=',NCOP,NSEG_SCR
                  call scramble(ixran,nseg_scr)
    c               print *,'IS,IXRAN(IS),ISHIFT=',is,ixran(is),ishift
    c             Last (shorter) segment
    c         Pad with random states, p(on)=percent on
                call randpx(1,ran)
    c     print *,'LASTFRAMEUSE,IFR=',lastframeuse,ifr
    c           if (auc(ifrr) .gt. nauc(ifrr)) write (6,7923)
    c    -        ifr,ifrr,ntodo,auc(ifrr),nauc(ifrr),itrack(ifr+ifrr)
    c7923       format(' ifr,ifrr=',2i6,' ntodo=',i6,' auc=',f12.1,
    c    -    ' nauc=',i6,' itrack=',i9)
    c     Normalize auc
    c      write (0006,9858) nframe2,nframeuse,lastframeuse
    c9858  format(' NFRAME2=',i6,' NFRAMEUSE=',i6,' LASTFRAMEUSE=',I6)
    c       if (auc(ifr) .lt. 0.0 .or. auc(ifr) .gt. 1.0) then
    c         print *,'IFR=',ifr,' AUC=',auc(ifr),' NAUC=',nauc(ifr),
    c    -     ' NFRAME2=',nframe2
    c         stop
    c       end if
            call trnsfr(auc_all(1,itr_ix),auc,nframe2)
          character*(*) lab,bondname(ntracks)
          common /aucw/ auc_all(MAXFRAMES21,MAXCOPY)
          call openfile(91,0,' ',1,'new',lab,llab,notfnd,0,1,1,0,0)
          close (91)
          character*(*) lab
          common /nnwork/ rmsd2d(MAX2D,MAX2D),fill(IFILL2)
          character* 132 line(maxrec)
    c     print *,'HBBRIDGEPRINT maxbridgelen,maxbridgemem=',
    c    -  maxbridgelen,maxbridgemem
          call zeroiti(lensum,0,maxbridgelen)
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
          character*1 ans
    c     print *,'BONDCORRSUM nbondcorr,nframe,nresbondcorr=',
    c    -  nbondcorr,nframe,nresbondcorr
          call trajlimtest(nframe,MAXFRAMES)
          call quiz(ans,icorrtyp,'b',' ',0,'bond correlation type',21,
          call quiz(ans,icorrtrans,'-',' ',0,
          call getint('Exponent of correlations',24,1,0,9,icorrexp,3)
          correxp=icorrexp
          call zeroiti(nng,0,MAX2D)
          call zeroiti(ifirstframe,0,MAX2D)
          call zeroiti(ing,0,MAX2D*MAX2D)
    c     if (iaggregate .gt. 0) then
    c       ihbtoresmax=0
    c       do ib=1,nbondcorr
    c         if (ihbtores(ib) .gt. ihbtoresmax) ihbtoresmax=ihbtores(ib)
    c       end do
    c     end if
            call readbitc(ires(1,ifr),ibframe,nbondcorr,30,MAXITEMS)
    c       Accumulate scalar product
    c                 Both on and off states are correlated
    c     Prepare the distance measure matrix
    c           Average states
    c          if (iaggregate .eq. 1)
    c     -      write (iout,8781) ib1,ib2,ing(ib1,ib2),ing(ib2,ib1),scp
    c8781      format(' IB1,2=',2i6,' ing(ib1,ib2)=',i6,
    c     -      ' ing(ib2,ib1)=',i6,' scp=',f10.5)
    c      write (iout,8791) (nng(i),i=1,ncorr)
    c      write (iout,8792) (ifirstframe(i),i=1,ncorr)
    c8791  format(' NNG:',10i6)
    c8792  format(' IFF:',10i6)
          character*(*) bondname
          character*132 line(maxrec)
          character*8 resnames(maxrec)
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
          common /bondpairs/ ihbpair(2,MAXBONDS),ihb_pair_res(3,MAXBONDS)
          character*1 corb
          character*8 aname(2),rname(2)
    c     print *,'BONDCORRP nclust,nhbcorr,iclust,bondname=',
    c    -  nclust,nhbcorr,iclust,bondname(1:lbondname)
          call write_traj_lim(iout,
          call zeroiti(ihistogram,0,10)
          corb=':'
    c      write (iout,9877) (indexbond(i),i=1,nhbcorr)
    c9877  format(' INDEXBOND:',/,(20i3))
              call findbestcorrep(iout,ifhbclust(ic),ilhbclust(ic),
                call lastchar(resnames(ixr1),lc1,8)
                call lastchar(resnames(ixr2),lc2,8)
          character*(*) yclab(nbonds)
          character*80 title
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
          call plotmat(ips,kc,rmsd2d,dc,nbonds,nbonds,0,0,0,0,1,1,
            call pswrite(ips,ixdel0,iytop,'m',1)
          call colcodeminmax(ips,25+ixcent,-iydel-5,0,ncolcode,
          character*(*) bondtype
          call zeroiti(nhbhist,0,10)
          character*8 atnames(maxrec),resnames(maxrsd)
          character*80 title,label2d(mx2d),label,remark
          character*(*) xtrajlab,inpfile
          character* 132 line(maxrec)
          character*(*) bondtype
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
          character*37 yclab(MAXBONDS)
          character*80 plotlab
          character*200 aucfile
    c     print *,'FINALIZEBOND nbondavg=',nbondavg
    c     Turn on sorting in clique clustering
          call zeroiti(nhbpers,0,MAXBONDS)
          call zeroiti(itrackf,0,MAXBONDS)
          call zeroiti(itrackl,0,MAXBONDS)
    c     print *,'FINALIZEBOND nbfound,nframe,ianchor2,nbresfound=',
    c    -  nbfound,nframe,ianchor2,nbresfound
    c     print *,'Number of ',bondtype(1:lbondtype),' bonds found=',nbfound
          call askyn('Do you want to end the track statistics at '//
          call trackstat(nbfound,nhbdist,nhbpers,maxlenon,maxlenoff,
          call printbonddist(nframe,nbfound,nhbdist,index2d,
          call filterbonds(n,nbfound,nhbdist,rhbdist,nhbpers,itrackf,
    c     cv contains the % of frames an atom is forming a bond
          call bondsum(nbfoundorig,nbfound,index2d,ixres,atnames,resnames,
          call makebondlab(1,nbfound,0,1,yclab,lyclab,irrix,ixresno,
            call getname(aucfile,laucfile,
            call print_auc(aucfile,laucfile,nbfound,nframe,yclab,lyclab)
          call plotbond(iw1,nbfoundorig,nbfound,nbresfound,nhbcorrclust,
    c     Plot the total # of bonds for successive frames
          call roundlimint(naabondmax+1,iyd2,nyd2)
          call roundlimint(nbondmax+1,iyd,nyd)
    c       write (iw0,*) i,' scr1=',scres(1,i),' dsum=',dsum
          call batchmean(nframe,0,data,'Average # of bonds',18,iw0,0,av,sd,
            call trnsfr(xtrajmod,xtraj,nframe)
          call plot2fun(iw1,nplot,xtrajmod,scres,scres,nfr,0.0,0.0,0,0.0,
    c       Calculate, plot and print the dot product of bonds
            call bondcorrsum(nbfoundorig,0,scpmin,scpmax,it1,it2,it3,
            call bondcorrprint(nbfound,iw0,line,index,iresno,index2d,
            call plotbondcorr(iw1,nbfound,yclab,lyclab,0,index2d,
            call clusterdistr(nbfound,iw0,rmsdlim,scpmin,scpmax,nhbdist,it1,
    c       print *,'Plotting the clustered ',bondtype(1:lbondtype),' bonds'
            call plotbond(iw1,nbfoundorig,nbfound,nbresfound,nhbcorrclust,
    c       print *,'Printing the clustered ',bondtype(1:lbondtype),' bonds'
            call bondcorrprint(nbfound,iw0,line,index,iresno,index2d,
            call askyn(
          call mapbondstorespairs(nbfound,nbfoundorig,nbresfilt,nbresfound,
    c     if (iauc+iresbondcorr .gt. 0)
    c    -  call condensetracks(nbfoundorig,nbresfilt,it1,it2,iw0,mxbonds)
          call condensetracks(nbfoundorig,nbresfilt,it1,it2,iw0,mxbonds)
          call res_res_bond(nres2d,nbfound,index2d,nhbdist,iresno,ixres,
          call plotbond(iw1,nbfoundorig,nbfound,nbresfilt,nhbcorrclust,
            call bondcorrsum(nbfoundorig,nbresfilt,scpmin,scpmax,it1,it2,
            call indexit(itemp1,1,nbresfilt,0)
            call makebondlab(1,nbresfilt,0,2,yclab,lyclab,irrix,ixresno,
            call plotbondcorr(iw1,nbresfilt,yclab,lyclab,1,itemp1,title,
            call clusterdistr(nbresfilt,iw0,rmsdlim,scpmin,scpmax,nhbdist,
            call bondcorrprint(nbresfilt,iw0,line,index,iresno,irrix,
            call plotbond(iw1,nbfoundorig,nbfound,nbresfilt,nhbcorrclust,
          call askyn(
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
    c     Condense the bond-bits to res-res bond bits in ires(*,nframe)
            call zeroiti(ibframe_rr,0,nbres)
            call readbitc(ires(1,ifr),ibframe,nbfoundorig,30,MAXITEMS)
            call savebitc(ires(1,ifr),ibframe_rr,nbres,30,MAXITEMS)
          character* 132 line(maxrec)
          character*(*) bondname
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
          common /bondpairs/ ihbpair(2,MAXBONDS),ihb_pair_res(3,MAXBONDS)
    c     print *,'MAPBONDSTORESPAIRS NBFOUND,NBFOUNDORIG,NBRESORIG=',
    c    -  nbfound,nbfoundorig,nbresorig
          call indexit(indexbond,1,nbfoundorig,0)
          call askyn(
            call getfiltlims(percmind,percmaxd,minresdistd,maxresdistd,
          call zeroiti(ihbtores,0,nbfoundorig)
          call zeroiti(nusepair,0,nbresorig)
          call indexit(itemp1,1,nbresorig,0)
          call indexit(indexres,1,nbresfound,0)
          call mrgsrti(6,indexres,itemp2,nbresfound,ifa,ila,itemp1,itemp3,
    c      write (iout,8711) (ihbtores(i),i=1,nbuse)
    c8711  format(' IHBTORES:',/,(15i5))
    c      write (6,8712) (irrix(i),i=1,nbresfound)
    c8712  format(' IRRIX:',/,(15i5))
          character*8 resnames(mxrsd)
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
          common /bondpairs/ ihbpair(2,MAXBONDS),ihb_pair_res(3,MAXBONDS)
    c     print *,'AUTOCORR_RES_RES NRESCOL=',nrescol,
    c    -  ' MXRSD,MXREC=',mxrsd,mxrec
          call auc_params(iauctype,lastframeinp,loffmin,nreusemax,
          call trackstat(nbres,nhbdist,nhbpers,maxlenon,maxlenoff,itrackf,
    c       Get the irr-th track
    c       write (iout,2000) irr,
    c    -    line(index(ia1))(irescol1:irescol2),iresno(ia1),isegno(ia1),
    c    -    line(index(ia2))(irescol1:irescol2),iresno(ia2),isegno(ia2),perc
            call persistence(nhbdist(irr),nhbpers(irr),itrackf(irr),
            call getbondtrack(irr,itrack,ifirstframe,lastframe,30,nframe)
            call autocorr(irr,irr,itrack,ifirstframe,lastframe,iframeunit,
          character*(*) trackfile
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
          common /bondpairs/ ihbpair(2,MAXBONDS),ihb_pair_res(3,MAXBONDS)
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
          character*200 trajnam,trajnam2
          common /trajname/ trajnam,trajnam2,ltrajnam,ltrajnam2,ifirsttraj,
            call zeroiti(itrack,0,nframe)
          close (iout_track)
          character*(*) trackfile
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
          common /bondpairs/ ihbpair(2,MAXBONDS),ihb_pair_res(3,MAXBONDS)
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
          character*200 trajnam,trajnam2
          common /trajname/ trajnam,trajnam2,ltrajnam,ltrajnam2,ifirsttraj,
          call blankout(trajnam,1,200)
          call lastchar(trajnam,ltrajnam,200)
            call zeroiti(itrack,0,nframe)
          call write_traj_lim(iout,' ',0,1,incr_tr,0)
          close (iout_track)
          character* 132 line(maxrec)
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
          common /bondpairs/ ihbpair(2,MAXBONDS),ihb_pair_res(3,MAXBONDS)
          corrmax=0.0
          corrmin=1.0
            call askyn('Do you want to connect correlated bonds',39,1,-1,
              call getreal(
              call askyn('Do you want to connect anti-correlated bonds',44,
                call getreal(
              call findCA(line,index,ifres(irr1),ilres(irr1),inamcol1,
              call findCA(line,index,ifres(irr2),ilres(irr2),inamcol1,
                cHe(k)=(c(k,iCA1)+c(k,iCA2))/2.0
    c          write (6,9782) nrr,irr1,irr2,iCA1,iCA2,nusepair(nrr)
    c9782      format(i4,' irr1,2=',2i5,' iCA1,2=',2i6,' nusepair=',i2)
    c       Connect bonds where the correlation distance measure is > corrmax
    c       or < corrmax
            call findCA(line,index,ifres(irr1),ilres(irr1),inamcol1,
            call findCA(line,index,ifres(irr2),ilres(irr2),inamcol1,
          character* 132 line(maxrec)
          character*8 aname
          common /nnwork/ rmsd2d(MAX2D,MAX2D),fill(IFILL2)
          character*(*) line
          common /askhex/ iaskhex(4),ishex(4)
          character*1 tab,ctrlM
          common /tab/ tab,ctrlM
          common /logging/ logfile,ipredict
          character*14 inplab(4)
    c     ihextyp=1: atom #; ihextyp=2: residue #; inphextyp=4: everything else
    c       Field was blank - return zero
    c       Interpret string as hexadecimal
    c       Check if hexadecimal
    c     Print the potential at the grid points
          common /nnwork/ phimap(MAXPHI,MAXPHI,MAXPHI)
          character*1 xyz
          common /axislab/ xyz(3)
              call getreal('X coordinate (999 to quit)',26,999999.0,x,0,0)
                call getreal('Y coordinate',12,999999.0,y,0,0)
                call getreal('Z coordinate',12,999999.0,z,0,0)
                call interpolate(x,y,z,gx,gy,gz,xstart,ystart,zstart,phi)
    c       Establish grid range to print
            call getreal('Grid X coordinate maximum',25,xstart+ngx*gx,xgmax,
            call getreal('Grid Y coordinate maximum',25,ystart+ngy*gy,ygmax,
            call getreal('Grid Z coordinate maximum',25,zstart+ngz*gz,zgmax,
    c           Mark grids to be dropped
    c                  write (77,7878) ia,ix,iy,iz,x,y,z,(c(k,ia),k=1,3),
    c     -              dd,rnear2,ndrop
    c7878              format(' ia=',i5,2x,3i5,' xyz=',3f10.4,' c=',3f10.4,
    c     -              ' dd,rnear2=',2f10.2,' ndrop=',i9)
              call fillinterpolate(phimap,ixgmin,ixgmax,iygmin,iygmax,
    c     print *,'FI nx1,nx,ny1,ny,nz1,nz=',nx1,nx,ny1,ny,nz1,nz
    c                 Empty line - stop for now
    c                 Empty line - stop for now
    c                 Empty line - stop for now
          common /nnwork/ phimap(MAXPHI,MAXPHI,MAXPHI)
            call interpolate(c(1,ia),c(2,ia),c(3,ia),gx,gy,gz,
            cv(ia)=phi
          common /nnwork/ phimap(MAXPHI,MAXPHI,MAXPHI)
          character*60 toplbl !ascii header
    c     print *,'READING MAP ',toplbl
    c     print *,'ivary, nbyte, intdat=',ivary, nbyte, intdat
          common /nnwork/ phimap(MAXPHI,MAXPHI,MAXPHI)
          character*(*) line
          character*(*) q,defname,name
          common /logging/ logfile,ipredict
          character*80 qline,ansline,defn
          call blankout(name,1,maxlen)
          call blankout(ansline,1,80)
          call lastchar(ansline,lc,80)
            call explanation(ihelp,0)
            call explanation(0,itip)
            call readquotestring(ansline,'"',i1,i2,ifail,80)
            call readquotestring(ansline,"'",i1,i2,ifail,80)
            call nextchar(ansline,ii,80)
            call nextblank(ansline,ii,80)
    c         Use default
          character*(*) line
          character*1 delim
          character*(*) line
          character* 132 line
          character*(*) q
          character*132 ansline,pline
          common /logging/ logfile,ipredict
    c     print *,'GETINT idef,noneg,maxval=',idef,noneg,maxval
    c     Allow idef == in
    c         Put default on query
          call blankout(ansline,1,132)
          call nextchar(ansline,ii,132)
              call explanation(ihelp,0)
          call nextblank(ansline,ii,132)
            call askyn('Do you want to interpret it as hexadecimal',42,1,-1,
          character*(*) q
          character*132 ansline,pline
          common /logging/ logfile,ipredict
          call blankout(ansline,1,132)
          call nextchar(ansline,ii,132)
              call explanation(ihelp,0)
          call lastchar(ansline,ilc,132)
          character*(*) q
          character*132 ansline,pline
          common /logging/ logfile,ipredict
    c       Put default on query
          call blankout(ansline,1,132)
          call nextchar(ansline,ii,132)
              call explanation(ihelp,0)
          call nextblank(ansline,ii,132)
          character*(*) line
          character*4 name
          call nextchar(line,ic,len)
          call nextblank(line,ic,len)
    c       Use trailing part of the name
          character*(*) q
          character*132 pline
          character*1 ans
          character*5 defans
          common /graphics/ npixx,npixy,maxpixx,maxpixy,idwmain,idwplot,
          common /logging/ logfile,ipredict
    c     idefans=-1: default no; idefans=+1: default yes
    c     iyn=1: yes -> ians=1, no -> ians=0
    c     iyn=0: yes -> ians=0, no -> ians=1
    c     print *,'ASKYN ihelp,itip=',ihelp,itip
              call explanation(ihelp,0)
              call explanation(0,itip)
            call lastchar(ans,ilc,1)
          character*(*) prefix,question0
          character*60 ansline
          character*80 line,line1(100),question
          character*1 char,defchar,ans,ch(2,100)
          character*60 q(100),promptlist,prompttype
          common /quizinfo/ nqst(600),lqst(600),iqfst(600),iqlst(600),
          common /logging/ logfile,ipredict
    c     if lprefix > 0 then the question will be composed of the prefix+question0
    c     if lprefix < 0 then the question will be composed of the question0+prefix
    c     drop the last nqincr questions
    c     print *,'QUIZ lprefix=',lprefix,' lqt0=',lqt0,' init=',init
    c       Initialize nq,lq,iqfst,iqlst
              call lastchar(promptlist(il),ilc,60)
    c       Quiz i has label prompptype(i) of length lqst(i)
    c       List of menu items in Quiz i: promptlist(iqfst(i))-promptlist(iqlst(i))
    c     Find the quiz label
    c       Combine prefix and question0
    c       Combine question0 and prefix
    c       Find signal character position (ich)
    c       Find last non-blank in the menu item line
            ch(1,iq)=q(iq)(ich:ich)
              call findcase(q(iq)(ich+2:ich+2),icase)
              call findcase(q(iq)(ich-2:ich-2),icase)
    c       Find lower-case of signal character
            call uplow(ch(1,iq),ch(2,iq),1,noabc)
    c       Generate <>-less menu item
              call uplow(line1(iq)(1:1),line1(iq)(1:1),2,noabc)
    c       Ask from the terminal - complete line with : and signal character
    c       Add help line
            ch(1,nqq)='?'
            ch(2,nqq)='?'
            call blankout(line,1,lq)
    c     Read and evaluate answer character
              char=ch(2,iq)
              char=ans
            call explanation(ihelp,0)
            call explanation(ihelp,0)
          character*60 promptlist,prompttype
          common /quizinfo/ nqst(600),lqst(600),iqfst(600),iqlst(600),
          character*55 submenu1(20)
          character*55 submenue(20)
          character*35 submenua(20)
          call findmenu('run type',8,ix1)
            call lastchar(submenu1(ix1-ix10),lc1,55)
              call findmenu(submenu1(ix1-ix10),lc1,ix2)
                  call lastchar(submenua(ix2-ix20),lc2,35)
                    call findmenu(submenua(ix2-ix20),lc2,ix3)
    c             Configuration edit submenus
                  call lastchar(submenue(ix2-ix20),lc2,55)
                    call findmenu(submenue(ix2-ix20),lc2,ix3)
          character*(*) label
          character*60 promptlist,prompttype
          common /quizinfo/ nqst(600),lqst(600),iqfst(600),iqlst(600),
          character*1 charin,charout
          character*1 abc,idig,digits,hexdigits
          common /charactersets/ ihex(25),abc(26,2),idig(10),digits(14),
    c     If iuptolow=1: up -> low; iuptolow=2: low -> up
              charout=abc(ic,3-iuptolow)
          character*1 charin
          character*1 abc,idig,digits,hexdigits
          common /charactersets/ ihex(25),abc(26,2),idig(10),digits(14),
    c     Returns icase=2: lower; 1: upper; 0: neither
          common /graphics/ npixx,npixy,maxpixx,maxpixy,idwmain,idwplot,
          common /rotmat/ matrot0(4,4),matrot(4,4),nomat0
          common /pbcrotmat/ torot_ac(3,3),torot_ca(3,3),tofac_ac,tofac_ca
          common /depthcuedat/ near,ifar,ramp0,idepth,idepthon,idrawh,
    c     print *,'DRAWPBC ioppbc=',ioppbc,' edge=',edgexyz,
    c    -  ' sizefac=',sizefac
    c     print *,'DRAWPBC cent=',cent
          call indexit(ixdup,1,48,0)
    c     Establish world coordinate system
    c       Rectangular cell
    c     Start polyhedron object definition
    c       Draw cell name
    c       print *,'wx,wxdr=',wx,wxdr,' wxlab,wylab=',wxlab,wylab
    c       Rectangular cell
    c       e2: half edges in pixel
    c       vp(k,1-4): y-z face at the -x side
    c       vp(k,5-8): y-z face at the +x side
            call shiftmol(vp,8,cent,vp,1.0)
            call writevertices(vp,8,ioppbc,ioutpdb,iconntyp,icrot,crot)
    c         Top polygon edge
              call connlinix(vp,m,mod(m,4)+1,ixdup,iconntyp,ioutpdb,48)
    c         Bottom polygon edge
              call connlinix(vp,m+4,mod(m,4)+5,ixdup,iconntyp,ioutpdb,48)
    c         Edge parallel to the prizm's axis
              call connlinix(vp,m,m+4,ixdup,iconntyp,ioutpdb,48)
    c       FCC
            call zeroit(vp,3*6)
            call shiftmol(vp,14,cent,vp,1.0)
            call finddup(vp,ixdup,14,nunique)
            call writevertices(vp,nunique,ioppbc,ioutpdb,iconntyp,
            call connlinix(vp,1,07,ixdup,iconntyp,ioutpdb,48)
            call connlinix(vp,1,08,ixdup,iconntyp,ioutpdb,48)
            call connlinix(vp,1,11,ixdup,iconntyp,ioutpdb,48)
            call connlinix(vp,1,12,ixdup,iconntyp,ioutpdb,48)
            call connlinix(vp,2,09,ixdup,iconntyp,ioutpdb,48)
            call connlinix(vp,2,10,ixdup,iconntyp,ioutpdb,48)
            call connlinix(vp,2,13,ixdup,iconntyp,ioutpdb,48)
            call connlinix(vp,2,14,ixdup,iconntyp,ioutpdb,48)
            call connlinix(vp,3,11,ixdup,iconntyp,ioutpdb,48)
            call connlinix(vp,3,12,ixdup,iconntyp,ioutpdb,48)
            call connlinix(vp,3,13,ixdup,iconntyp,ioutpdb,48)
            call connlinix(vp,3,14,ixdup,iconntyp,ioutpdb,48)
            call connlinix(vp,4,07,ixdup,iconntyp,ioutpdb,48)
            call connlinix(vp,4,08,ixdup,iconntyp,ioutpdb,48)
            call connlinix(vp,4,09,ixdup,iconntyp,ioutpdb,48)
            call connlinix(vp,4,10,ixdup,iconntyp,ioutpdb,48)
            call connlinix(vp,5,07,ixdup,iconntyp,ioutpdb,48)
            call connlinix(vp,5,09,ixdup,iconntyp,ioutpdb,48)
            call connlinix(vp,5,13,ixdup,iconntyp,ioutpdb,48)
            call connlinix(vp,5,11,ixdup,iconntyp,ioutpdb,48)
            call connlinix(vp,6,08,ixdup,iconntyp,ioutpdb,48)
            call connlinix(vp,6,10,ixdup,iconntyp,ioutpdb,48)
            call connlinix(vp,6,14,ixdup,iconntyp,ioutpdb,48)
            call connlinix(vp,6,12,ixdup,iconntyp,ioutpdb,48)
    c       Hexagonal prism
            call shiftmol(vp,12,cent,vp,1.0)
            call writevertices(vp,12,ioppbc,ioutpdb,iconntyp,icrot,crot)
    c         Top polygon edge
              call connlinix(vp,m,mod(m,6)+1,ixdup,iconntyp,ioutpdb,48)
    c         Bottom polygon edge
              call connlinix(vp,m+6,mod(m,6)+7,ixdup,iconntyp,ioutpdb,48)
    c         Edge parallel to the prizm's axis
              call connlinix(vp,m,m+6,ixdup,iconntyp,ioutpdb,48)
    c       Truncated octahedron cell
            call zeroit(vp,(3*6*8))
    c       Hexagon face (+x,+y,+z) quadrant
    c       Hexagon face 2 (+x,+y,-z) quadrant
    c       Hexagon face 3 (-x,+y,-z) quadrant
    c       Hexagon face 4 (-x,+y,+z) quadrant
    c       Hexagon face 5 (+x,-y,+z) quadrant
    c       Hexagon face 6 (+x,-y,-z) quadrant
    c       Hexagon face 7 (-x,-y,-z) quadrant
    c       Hexagon face 8 (-x,-y,+z) quadrant
            call shiftmol(vp,48,cent,vp,1.0)
            call finddup(vp,ixdup,48,nunique)
            call writevertices(vp,nunique,ioppbc,ioutpdb,iconntyp,
                call connlinix(vp,m+(ifa-1)*6,mod(m,6)+1+(ifa-1)*6,ixdup,
    c       print *,'t,hh,hl=',t,hh,hl
    c       x-y coordinates of wrapping hexagon vertices
            call zeroit(vp,6)
    c       Top rhombuses
    c       First the chair hexagon
              call trnsfr(vphu(1,i),hxy(1,i),2)
            call trnsfr(vp(1,3),vphu,18)
    c       Bottom rhombuses
            call trnsfr(vphl(1,1),vphu(1,1),21)
            call trnsfr(vp(1,9),vphl,18)
            call finddup(vp,ixdup,14,nunique)
            call writevertices(vp,nunique,ioppbc,ioutpdb,iconntyp,
    c       Top rhombuses
              call connlinix(vp,m+2,mod(m,6)+3,ixdup,iconntyp,ioutpdb,48)
              call connlinix(vp,1,k+2,ixdup,iconntyp,ioutpdb,48)
    c       Bottom rhombuses
              call connlinix(vp,8+m,mod(m,6)+9,ixdup,iconntyp,ioutpdb,48)
              call connlinix(vp,2,k+8,ixdup,iconntyp,ioutpdb,48)
    c       Side connections
              call connlinix(vp,i+2,i+8,ixdup,iconntyp,ioutpdb,48)
    c       Octahedral - draw parallelepiped
    c       print *,'DRAWPBC ioppbc=9 '
            call arrsum(vp(1,1),edge_gen(1,1),vp(1,2),3)
            call arrsum(vp(1,2),edge_gen(1,2),vp(1,3),3)
            call arrsum(vp(1,1),edge_gen(1,2),vp(1,4),3)
            call arrsum(vp(1,1),edge_gen(1,3),vp(1,5),3)
            call arrsum(vp(1,5),edge_gen(1,1),vp(1,6),3)
            call arrsum(vp(1,6),edge_gen(1,2),vp(1,7),3)
            call arrsum(vp(1,5),edge_gen(1,2),vp(1,8),3)
            call writevertices(vp,8,ioppbc,ioutpdb,iconntyp,
              call connlinix(vp,m,mod(m,4)+1,ixdup,iconntyp,ioutpdb,48)
              call connlinix(vp,m+4,mod(m,4)+5,ixdup,iconntyp,ioutpdb,48)
              call connlinix(vp,m,m+4,ixdup,iconntyp,ioutpdb,48)
    c        write (6,7171) (i,(vp(k,i),k=1,3),i=1,nv)
    c7171    format(i5,3f10.3)
          character*3 pbcres
          common /pbcresname/ pbcres(10)
    c     print *,'WRITEVERTICES n=',n,' icrot=',icrot
    c     print *,'crot=',crot
              c8(k)=vp(k,ia)
            call writepdbd(ioutpdb,c8,ia,1,'VERT',pbcres(ioppbc),'V',1.0,0.)
    c     Just draw a line connecting p(i1) & p(i2)
    c     iconntyp = 1,2,3: SGI, PDB, both
    c         Duplicate - delete
              call trnsfr(p(1,i-ndel),p(1,i),3)
    c      write (6,7777) n,nunique,(ixdup(i),i=1,n)
    c7777  format(' FINDDUP n,nunique=',2i6,(' ixdup=',10i3,3x,10i3))
          common /depthcuedat/ near,ifar,ramp0,idepth,idepthon,idrawh,
    c     print *,'SETCOLOR idepthon,ic,idepth=',idepthon,ic,idepth
          common /depthcuedat/ near,ifar,ramp0,idepth,idepthon,idrawh,
    c     ref. Computers in Chemistry Vol 13, No 3, pg 191, 1989
    c     Approach: Construct a vector A from ca atom i to i-1. Construct B
    c     from i to i+1.  Find V1, the vector which bisects A and B.  V1 is
    c     perpendicuar to the helix axis (for a perfect mathematical helix).
    c     Let i = i + 1 and repeat the procedure, finding V2.  As V1 and V2
    c     are both perp. to the helical axis, their cross product gives the
    c     helix direction.   Average over all possible tetrads of ca atoms,
    c     or use a 3-D fitting method to give the direction based on the
    c     points calculated to lie on the axis.
    c
    c     P1 and P2 are the vectors from the origin to CA 1 and CA 2.
    c
    c     The radius is a calculated as:
    c     |dH|**2 - |p2-p1|**2
    c     r = --------------------
    c     2 * |(p1-p2) dot v2|
    c     where d= (p2-p1) dot h
    c
    c     H1 and H2 are the position vectors
    c     H1 = P1 + r*V1
    c     H2 = P2 + r*V2
    c
          character*60 Message
    c     Arrays of length MAXHX
    c     Arrays of length 3,2*MAXHX
    c     load p1 and p2
             call dvset(p1,co(1,i))
             call dvset(p2,co(1,i+1))
    c     get vector v1
             call dvdif (p1,co(1,i-1),a)
             call dvdif (p1,co(1,i+1),b)
             call dvnorm(a)
             call dvnorm(b)
             call dvsum(a,b,v1)
             call dvnorm(v1)
    c     get vector v2
             call dvdif (p2,co(1,i),a)
             call dvdif (p2,co(1,i+2),b)
             call dvnorm(a)
             call dvnorm(b)
             call dvsum(a,b,v2)
             call dvnorm(v2)
    c     H=direction of axis
             call dcross(v1,v2,h)
             call dvnorm(h)
             call dvsum(h,hsum,hsum)
    c     calculate radius
             call dvdif (p1,p2,p1mp2)
             call dvdif (p2,p1,p2mp1)
             call dvmul(h,d,tmp)
    c     kahn
    c     calculate the points on the axis
             call dvmul(V1,r,tmp)
             call dvsum(P1,tmp,H1(1,hcount))
             call dvmul(V2,r,tmp)
             call dvsum(P2,tmp,H1(1,hcount+1))
    c     hcount=hcount+1
          call dvset(dir,hsum)
          call dvnorm(dir)
          call parlsq(h1,hcount-1,docircfit,dir,ip,fp,rms,0)
             call circfit(co,nats,dir,ip)
    c     adjust ip to be next to the first alpha carbon, not the second
             call dvdif (co(1,1),ip,tmp)
             call dvproj(dir,tmp,tmp)
             call dvsum(ip,tmp,ip)
          call RMScalc(co,nats,dir,ip,fp,RMS)
          call writeout_h(dir,ip,fp,rms,message,iprint)
    c     from kahn's package
    c     from kahn's package
    c     calculates the axis of a helix from alpha carbon coordinates using
    c     linear least squares regressions of atom number versus x,y, and z
    c     coordinates.
          character*60 message
    c     print *,'PARLSQ MAXHX,n=',MAXHX,n
    c      write (6,7777) ((co(k,i),k=1,3),i=1,n)
    c7777  format(' PARLSQ input:',/,(3f10.4))
    c      print *,'docircfit=',docircfit
    c     use circular fit for ip
             call circfit(co,n,dir,ip)
    c     use least-squares initial point
          call dvnorm(dir)
          call RMScalc(co,n,dir,ip,fp,RMS)
          call writeout_h(dir,ip,fp,rms,message,iprint)
    c     Modified from the original to calculate the sd of the atom to axis dist
    c     calculate rms deviations
             call dvdif (co(1,i),ip,tmp)
             call dvproj(dir,tmp,tmp)
             call dvsum(ip,tmp,fp)
          character*60 message
    c     don't waste time
    c     Array of 3,2*MAXHX length
    c     calculates the best-fit circle to a set of data and returns the
    c     center in ip
    c     x=the coordinates
    c     dir(input) the direction vector
    c     ip(output) the position vector
    c     print *,'CIRCFIT nats,MAXHX=',nats,MAXHX
    c     rotate the coordinates
    c
    c     read input
          call dvset(ax(1),ip)
          call dvset(ax(4),dir)
    c      write(*,*) 'got ax:'
    c      write(*,'(6(x,f8.3))') ax
    c      write(*,*) 'got nats:',nats
    c      write(*,*) 'got coords:'
    c      do i=1,nats
    c         write(*,'(3(x,f8.3))') (x(k,i),k=1,3)
    c      end do
          call polar(dir(1),dir(2),dir(3),r,thet,fi)
    c     write(*,*) 'Theta,phi',thet,fi
             call dvset(xp(1,iat),x(1,iat))
             call rotabout(xp(1,iat),zero,-thet,'z')
             call rotabout(xp(1,iat),zero,-fi,'y')
          call dvset(tmp,dir)
    c     write(*,'(''sum'',5(x,f8.3))') sum
          call rotabout(xp0,zero,fi,'y')
          call rotabout(xp0,zero,thet,'z')
          call dvset(ax,xp0)
    c     save output
          call dvset(ip,ax)
    c     write(*,*) 'Output ip, dir'
    c     write(*,'(6(x,f8.3))') ax
    c     Return polar coordinates for point (x,y,z)
    c     NOTE: Unlike most conventions, phi is the angle between R and the Z axis
             cosa=dble(z)/dble(r)
    c     rotate vector v about point v0 by t radians along the axis
    c     specified.
          character*1 axis
          ct=cos(t)
    c     /*- transpose an NxN matrix rot -*/
    c     /*- z =  x / y -*/
    c     /*- sets vector a=b -*/
    c     /*- normalize the vector x -*/
    c        write(*,*) 'vnorm: can''t normalize zero vector'
             continue
    c     /*- normalize the vector x -*/
    c     /*- compute z = x cross y -*/
    c     /*- compute z = x - y -*/
    c     /*- z = x * y -*/
    c     /*- z =  x / y -*/
    c     /*- z = x + y -*/
    c     /*- z = x + y -*/
    c/*-
    c     compute z = the projection on x of y
    c     projxy requires that x be a unit vector
    c     formula: projxy = (y dot x)*x
    c-*/
          call dvset(tmp,x)
          call dvnorm(tmp)
          call dvmul(tmp,ddot(tmp,y),z)
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
    c     Calculate the relaxation time of the ires-th stored residue
    c     print *,'RESAUTOCORR ires,nframe,increment=',ires,nframe,increment
          call zeroit(xyplot,2*ncf)
          character*8 resnames
          character*3 hphres
          character*1 ans
    c     EIS KD WHI
            call quiz(ans,ihydtyp,' ',' ',0,
                call getreal('Residue '//hphres(i)//' hydrophobicity',26,
            cv(ia)=hph3(ihydtyp,ir)
    c     Set solvent atoms to unknown
            cv(ia)=hph3(ihydtyp,23)
    c     Calculating the circular variance for solute atoms and solvent molecules
    c     print *,'CVLIST n,nslt,naslv=',n,nslt,naslv
          call zeroit(rprox,nslt)
              cv(ia)=1.d0-dsqrt(xnum*xnum+ynum*ynum+znum*znum)/den
              cv(ia)=0.0
    c       CV for solvents
              cvi=0.0
                cv(i)=cvi
              call trnsfr(ctemp(1,nslt+1),c(1,nslt+1),3*(n-nslt))
              call mrgsrt(6,ixtemp,rtemp,numslv,ifa,ila,it1,ctemp,n)
                  cv(i)=rtemp(ia)
                  call trnsfr(c(1,i),ctemp(1,ix),3)
          character*8 repnam
          character*(*) title
          common /nnwork/ mx(MAXCV,MAXCV),rij(3,MAXCV,MAXCV),
          character*4 yclab(1)
          character*47 maptyp(2)
          character*8 atnam
          character* 132 line(maxrec)
          common /graphics/ npixx,npixy,maxpixx,maxpixy,idwmain,idwplot,
          call indexit(ixshuffle,1,MAXCV,0)
    c     print *,'inamcol1,inamcol2,lnam=',inamcol1,inamcol2,lnam
    c     print *,'renam=',repnam(1:lnam),'*'
            call leftadjustn(atnam,atnam,8)
    c       write (77,*) ia,iresn,' atnam=',atnam(1:lnam),'*'
    c         New residue
              call trnsfr(ca(1,nres+1),c(1,ia),3)
    c           Normalize the rij vectors
    c     Internally 1-CV is stored/calculated
          cvavfor(1)=0.0
          cvavfor(nres)=0.0
          cvavback(1)=0.0
          cvavback(nres)=0.0
            call zeroitd(rijsum,3)
              cv=dsqrt(rijsum(1)**2+rijsum(2)**2+rijsum(3)**2)/dijsum
              cvrow(iy)=cv
    c     Now scan backwards
            call zeroitd(rijsum,3)
              cv=dsqrt(rijsum(1)**2+rijsum(2)**2+rijsum(3)**2)/dijsum
              cvcol(iy)=cv
            cvavfor(ix)=forsum(ix)/(nres-ix)
            cvavback(ix)=backsum(ix)/(ix-1)
            cvav(ix)=(forsum(ix)+backsum(ix))/nres
    c      do i=1,10
    c        write (6,7733) i,(mx(i,j),j=1,10)
    c      end do
    c7733  format(' MX i=',i3,' mx(i)=',10i3)
    c     Plot map
          call plotmat(ips,mx,rmx,dmx,nres,nres,0,0,0,0,1,nrep,ixdel,iydel,
    c     Plot row averages for domain map
          call colstrip(cvavfor,nres,ixdel,scalefac,iydown,nrep,' F',
          call colstrip(cvavback,nres,ixdel,scalefac,iydown+18,nrep,' B',
          call colstrip(cvav,nres,ixdel,scalefac,iydown+36,nrep,' T',
          call colcode01(ips,ixdel,iydown+5,ncolcode,nrep)
    c     Plot the type of the plot and close the plot/page
            call rgbcolor(ips,9)
            call psshow(ips,maptyp(imaptyp),47)
            close (iwr)
            close (ips)
          character*(*) title,title2,xylab,yclab(nyclab)
    c     kc,rc,dc: matrix in integer, real*4 and real*8 formats
    c     mxr,mxrr,mxrd: dimensions of kc,rc,dc (only one of them can be > 1)
    c     nresx, nresy: matrix dimension
    c     ishufflex, ixinc: X axis label of the ith row is ixshufflex(i+ixinc)
    c     ishuffley, iyinc: Y axis label of the ith row is ixshuffley(i+iyinc)
    c     ixincshift,iyincshift: deduct from the matrix indices
    c     navg: number of residues the matrix entries were averaged over
    c     nrep: When >1 then it was a repeat plot
    c     ixdel,iydel: increments in the x and y direction
    c     scalefac: Overall scaling factor (500*550 matrix fills the page
    c     with scalefac=1
    c     rcmin,rcmax: min and max of matrix value allowed
    c     ncol, maxcol: requested and maximum number of colors
    c     inc: increment per matrix row/col
    c     print *,'PLOTMAT ips,ncol=',ips,ncol,' nresx,nresy=',nresx,nresy
    c     print *,'PLOTMAT mxr,mxrr,mxrd=',mxr,mxrr,mxrd
    c     print *,'PLOTMAT mx,mxix=',mx,mxix
    c     print *,'PLOTMAT iydel,iydel0,idown,iyinc=',
    c    -  iydel,ixdel0,idown,iyinc
    c     print *,'PLOTMAT ixdel,ixdel0,iytop=',iytop,ixdel,ixdel0
    c     print *,'PLOTMAT ltitle,ltitle2=',ltitle,ltitle2,
    c    -  ' xylab=',xylab(1:lxylab)
    c      do i=1,10
    c        write (6,7733) i,(kc(i,j),j=1,10)
    c      end do
    c7733  format(' KC i=',i3,' kc(i)=',10i3)
    c     if (ltitle .gt. 0) print *,'PLOTMAT title=',title(1:ltitle)
    c     if (ltitle2 .gt. 0) print *,'PLOTMAT title=',title2(1:ltitle2)
    c     print *,'PLOTMAT ITRAJNAME=',itrajname
    c     write (6,8422) 'X',(ixshufflex(i),i=1,nresx)
    c     write (6,8422) 'Y',(ixshuffley(i),i=1,nresy)
    c8422 format(' IXSHUFFLE',a1,':',/,(20i4))
            call pswrite(ips,ixdel,iytop,'m',1)
            call psshow(ips,title,ltitle)
    c       Plot both names
            call pswrite(ips,ixdel,iytop,'m',1)
            call write_traj_lim(ips,' ',1,2,incr_tr2,1)
            call pswrite(ips,ixdel,iytop,'m',1)
            call write_traj_lim(ips,' ',1,1,incr_tr,1)
            call pswrite(ips,ixdel,iytop,'m',1)
            call write_traj_lim(ips,' ',1,itrajname,incr_tr,1)
            call pswrite(ips,ixdel,iytop,'m',1)
            call psshow(ips,title2,ltitle2)
    c     print *,'PLOTMAT iydel,ixinc,iyinc,inc=',
    c    -  iydel,ixinc,iyinc,inc
    c       write (77,*) 'iyy,iyinc,iy=',iyy,iyinc,iy
    c           Convert rc(ix,iy) or dc(ix,iy) into color code
    c          write (77,7688) ix,iy,rc(ix,iy),kccurr
    c7688      format('ix,iy=',2i6,' rc=',f10.5,' kcurr=',i3)
    c               Close the rectangle that is open
                    call pswrite(ips,ixdelmat+(ixx-1)*inc,iydel+(iyy)*inc,
                    call pswrite(ips,ixdelmat+(ixx-1)*inc,iydel+(iyy-1)*inc,
                    call pswrite(ips,ix0-inc,iy0-inc,'l',1)
    c               Open new rectangle
                      call rrgbcolor(ips,kccurr,100,0)
                      call rgbcolor(ips,kccurr)
    c                write (77,8766) ixx,iyy,ix0,iy0,inc
    c8766  format(' ixx,iyy=',2i4,' ix0,iy0=',2i4,' inc=',i2)
                    call pswrite(ips,ix0-inc,iy0-inc,'m',1)
                    call pswrite(ips,ix0-inc,iy0,'l',1)
    c         Close last opened rectangle
              call pswrite(ips,ixdelmat+(nresx)*inc,iydel+(iyy)*inc,'l',1)
              call pswrite(ips,ixdelmat+(nresx)*inc,iydel+(iyy-1)*inc,'l',1)
              call pswrite(ips,ix0-inc,iy0-inc,'l',1)
    c     Write residue/frame number scale
    c     print *,'Matrix entries plotted'
    c     Plot protein id
    c     if (nrep .le. 1) then
    c       call rgbcolor(ips,9)
    c       call pswrite(ips,ixdel,iydel-idown1,'m',1)
    c     end if
          call rgbcolor(ips,9)
            call pswrite(ips,ix,iydel-idown-iydown,'m',1)
            call roundlim(xmax,xdiv,nxdiv)
            call roundlim(ymax,ydiv,nydiv)
            call roundlimint(nresx,ixdiv,nxdiv)
            call roundlimint(nresy,iydiv,nydiv)
    c     print *,'DRAWRECT nresy,iydiv,iydel=',nresy,iydiv,iydel
    c     print *,'DRAWRECT ixdelmat,inc=',ixdelmat,inc
          call drawrect(ipsdraw,9,1,ixdelmat,ixdelmat+nresx*inc,iydel,
    c       write (6,8768) nresx,irinc,xdiv,xdiv_i,index(nresx)
    c8768   format(' nresx=',i6,' irinc=',i6,' xdiv,xdiv_i=',2f6.2,
    c    -    ' index(nresx)=',i6)
    c           Upper tick
                call pswrite(ips,ixdelmat+ir,iydel+nresy*inc+idown,'m',1)
                call pswrite(ips,ixdelmat+ir,iydel+nresy*inc+idown1,'l',1)
    c           Lower tick
                call pswrite(ips,ixdelmat+ir,iydel-idown,'m',1)
                call pswrite(ips,ixdelmat+ir,iydel-idown1,'l',1)
    c           Scale
                call pswrite(ips,ixdelmat+ir-ixsshift,iydel-idown2,'m',1)
    c    -               xdiv .ge. 100.0) then
    c         Left tick
              call pswrite(ips,ixdelmat,iydel+ir,'m',1)
              call pswrite(ips,ixdelmat-(idown1-idown),iydel+ir,'l',1)
    c         Right tick
              call pswrite(ips,ixdelmat+nresx*inc,iydel+ir,'m',1)
                call pswrite(ips,ixdelmat+nresx*inc+(idown1-idown),iydel+ir,
                call pswrite(ips,ixdelmat+nresx*inc-(idown1-idown),iydel+ir,
    c           Scale
                call pswrite(ips,ixdel0,iydel+ir-4,'m',1)
    c    -               ydiv .ge. 100.0) then
    c     print *,'PLOTMAT NYCLB,LYCLAB=',nyclab,lyclab
    c       Print residue identifiers
              call pswrite(ips,ixdelmat+nresx*inc,iydel+(iy-1)*inc+inc/3,
          character*9 filename
    c     print *,'CONTRACTMAT nxo,nyo,isort=',nxo,nyo,ndim,isort
          call getint(
    c     print *,'NAVG=',navg
    c         Sort columns first
                call trnsfr(rij(1,iy),temp,nxo)
              call trnsfi(itemp1,ixy,nyo)
                call trnsfr(temp,rij(1,itemp1(iy)),nxo)
                call trnsfr(rij(1,itemp1(iy)),rij(1,iy),nxo)
                call trnsfr(rij(1,iy),temp,nxo)
          character*(*) lab
    c     Writes a PS command with the minimum number of blanks
          character*80 line
          call writeint(line,icol,i1,lenw)
          call writeint(line,icol,i2,lenw)
          character*(*) line
    c     Writes the integer int into line starting at icol
    c     Draw a rectangle
    c      write (6,1000) kc,ix0,ix1,iy0,iy1
    c1000  format(' DRAWRECT kc=',i1,' ix1,2=',2i5,' iy1,2=',2i5)
            call rgbcolor(ips,kc)
            call pswrite(ips,ix0,iy0,'m',1)
            call pswrite(ips,ix0,iy1,'l',1)
            call pswrite(ips,ix1,iy1,'l',1)
            call pswrite(ips,ix1,iy0,'l',1)
            call pswrite(ips,ix0,iy0,'l',1)
          character*(*) label
            call rrgbcolor(iplot,i,100,1)
          call rrgbcolor(iplot,i,1,0)
    c       Use rcmin,rcmax
    c       Use 0 - n/rn
          call psshow(iplot,label,llabel)
          character*2 lab
              call rgbcolor(ips,kcprev)
            call rgbcolor(ips,kcprev)
            call rgbcolor(ips,9)
            call psshow(ips,lab,2)
    c     Plot color code
              call rgbcolor(ips,9)
              call rgbcolor(ips,icol)
            call drawrect(-ips,9,1,ix0+35,ix0+45,-iyd2,-iyd0,nrep)
          character*6 limline
    c     Plot color code
    c     print *,'COLCODEMINMAX ips,ncode,iydown=',ips,ncode,iydown
              call rgbcolor(ips,9)
              call rgbcolor(ips,icoldrw)
            call rgbcolor(ips,9)
          call rgbcolor(ips,9)
            call drawrect(-ips,9,1,ix0+35,ix0+45,-iyd2,-iyd0,nrep)
    c     Set the PS color to the SGI colors (0-7); reverse 0 & 7
    c     0:white; 1:red; 2:green; 3:yellow; 4:blue; 5:pink; 6; cyan; 7:orange
    c     8:yellow-green 9:back
    c     Set the PS color to the rainbow colors (0-7);
    c     0:white; 1:red; 2: orange; 3: yellow; 4:green;
    c     5:cyan; 6:blue; 7:violet; 8:indigo; 9:black 10: magenta
          common /colorinfo/ ncolcode,maxcolcode
    c     Set the color on a continous 0-1 scale
    c       rcol=2.0*float(icol-1)/float(imax-1)
    c       if (rcol .le. 1.0) then
    c         write (iout,1000) (1.0-rcol),rcol,0.0
    c       else
    c         rcol=rcol-1.0
    c         write (iout,1000) 0.0,(1.0-rcol),rcol
    c       end if
    c*****Postscript plot of the development of DSSP assignments
          character*(*) itit,xlab,ylab
          character*1 typc
          character*21 ssname
          common /dsspnames/ lssname(9),ssname(9),typc(9)
    c     nxd, nyd: increment between tics in the x, y axes, resp
    c     itit: array containing the title; ntit: number of characters in itit
    c     xlab, ylab: labels of the x, y axes
    c     nxlab, ylab: number of characters in the labels of the x,y axes
    c     print *,'PLOT nss,icx,maxconf=',nss,icx,maxconf
    c     print *,'PlOT nxd,nyd,nframe=',nxd,nyd,nframe
    c     xm0=0.075*xm
    c     ym0=0.09*ym
    c     print *,'ixd,iyd,xmin,xmax=',ixd,iyd,xmin,xmax
            call psheader(ips,itit,ntit,-30,-130,830,830,1,ipspage)
    c       write (ips,1001) xm0+0.3*xmm,ym0+1.03*ymm,' m'
            call psshow(ips,itit,ntit)
    c       write (ips,1001) xm0+0.45*xmm,ym0-0.10*ymm,' m'
            call psshow(ips,xlab,nxlab)
            call psshow(ips,ylab,nylab)
              call rgbcolor(ips,is)
    c         write (ips,1001) xm0+shiftl(is)*xmm,ym0-0.15*ymm,' m'
              call psshow(ips,ssname(is),lssname(is))
              call psshow(ips,':',1)
            call rgbcolor(ips,9)
                call rgbcolor(ips,itypss(iss))
          character*(*) title,bondname
          character* 132 line(maxat)
          character*200 trajnam,trajnam2
          common /trajname/ trajnam,trajnam2,ltrajnam,ltrajnam2,ifirsttraj,
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
          common /bondpairs/ ihbpair(2,MAXBONDS),ihb_pair_res(3,MAXBONDS)
          character*1 ans
          character*13 yclab(40)
          character*80 question
    c     print *,'RES_RES nframe,NBFOUND,MAXBONDS=',nframe,nbfound,MAXBONDS
          call indexit(ixshuffle,1,MAX2D,0)
          call zeroiti(irrix,0,nres)
          call zeroiti(it2,0,nres)
          call zeroiti(it3,0,nres)
    c       write (iout,9876) i,ihb,(ihbpair(k,ihb),k=1,2),
    c    -    (ixres(ihbpair(k,ihb)),k=1,2)
    c9876   format(' I,IHB=',2i5,' IHBPAIR=',2i6,' IXRES=',2i6)
    c     irrix(ir) is one if residue # ir is in a bond listed
    c     it2 points from the original list to the short list
    c     it3 points from the short residue list to the original list
          call askyn(question,lq,1,+1,ireduceplot,73,0)
          call quiz(ans,iresres,'c',' ',0,'treatment of contacts',21,
            call indexit(it2,1,nres,0)
            call indexit(it3,1,nres,0)
    c     Calculate cumulative res-res bond percentage
            call writesmall(iout,a1,nreshb,jr,it3(jr))
    c     Calculate average res-res bond percentage
            call zeroit(a1,nreshb)
            call writesmall(iout,a1,nreshb,jr,it3(jr))
    c     Generate contact statistics ignoring multiple contacts
    c     Sort contact list
          call sortbondlist(ixres,nbfound,indexbond,maxrsd,mxbonds)
            call readbitc(ires(1,ifr),ibond,nbfound,30,MAXITEMS)
            call zeroit (a1,nreshb)
            call writesmall(iout,a1,nreshb,jr,it3(jr))
    c       call indexit(irrix,1,nreshb,0)
            call indexit(irrix,1,n_res_res,0)
            call getreal(question,lq,rhbmax,rhbscalemax,0,74)
            call getint('Last residue to plot on the Y axis',34,
            call getint('First residue to plot on the Y axis',35,
    c         Gather residue names and numbers
    c           print *,iy,' YCLAB=',yclab(iy-iyinc)(1:lyclab)
            call plotmat(iplot,ing,rmsd2d,dc,nreshb,nreshby,0,0,0,iyinc,
    c       Draw boundary between solute molecules
            call rgbcolor(iplot,9)
    c           print *,'BOUND at i=',i,' ir=',it3(i),' ia=',ifres(it3(i+1))
            call colcodeminmax(iplot,25+ixcent,-iydel-5,0,ncolcode,
    c1015  format('( Upper-left triangle: iy is donor; ',
    c     -  'Lower-right triangle: ix is donor ) show')
          character*80 line
    c     print *,'WRITESMALL n,iy,iyorig=',n,iy,iyorig
            call blankout(line,1,80)
    c     Sort a list of positive integers
    c     print *,'MRGS n=',n
    c     print *,'LIST=',(list(i),i=1,n)
          call indexit(it1,1,n,0)
    c     do i=1,n
    c       t1(i)=list(i)
    c     end do
          call mrgsrti(6,it1,list,n,it2,it3,it4,it5,n)
    c     print *,'T=',(t1(i),i=1,n)
    c     do i=1,n
    c       list(i)=t1(i)
    c     end do
    c     print *,'LIST=',(list(i),i=1,n)
          character*3 lab
    c     Sort a list of distinct positive integers of maximum value maxval
    c     Use it for long lists; listlen about the same as maxval
          call zeroiti(itemp,0,maxval)
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),fill(IFILL6),
          common /bondpairs/ ihbpair(2,MAXBONDS),ihb_pair_res(3,MAXBONDS)
          call mrgsrt(6,indexbond,a1,listlen,i2,i3,i4,a2,listlen)
    c      do i=1,listlen
    c        ii=indexbond(i)
    c        iii=ixres(ihbpair(1,ii))*MAXBONDS+ixres(ihbpair(2,ii))
    c        write (6,5511) i,ii,(ixres(ihbpair(k,ii)),k=1,2),iii
    c5511    format(i4,' ixb=',i3,' ir1,2=',2i4,' a=',i9)
    c      end do
    c     Get the atom and residue bond sums
          character*8 atnames(maxrec),resnames(maxrsd)
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
          common /bondpairs/ ihbpair(2,MAXBONDS),ihb_pair_res(3,MAXBONDS)
    c     print *,'BONDSUM nhbtot,nhbfilt=',nhbtot,nhbfilt
          call askyn(
          call zeroiti(nhb_atot,0,maxat)
          call zeroiti(nhb_rtot,0,maxres)
            call readbitc(ires(1,ifr),ihb,nhbtot,30,MAXITEMS)
    c       Get the atom and residue bond sums
            call zeroiti(ihb_a,0,maxat)
            call zeroiti(ihb_r,0,maxres)
                call lastchar(atnames(ib1),lc,8)
                call lastchar(atnames(ib2),lc,8)
                call lastchar(resnames(ir1),lc,8)
                call lastchar(resnames(ir2),lc,8)
    c*****Plot the time-course of hydrogen, hydrophobic bonds & salt bridges
          character*(*) itit,xlab,ylab,bondname
          character*8 atnames(maxrec),resnames(maxrsd)
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
          character*41 clstyp
          common /cluster_typ/ nclstyp,inumclst(9),ireadcutoff(9),
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
          character*1 corrtrans(3)
          character*36 corrtyp(3)
          character*37 bondlab(1)
    c     nxd, nyd: number of tics in the x, y axes, resp
    c     itit: array containing the title; ntit: number of characters in itit
    c     xlab, ylab: labels of the x, y axes
    c     nxlab, ylab: number of characters in the labels of the x,y axes
    c     ibondlab=1: bond label; ibondlab=2: res pair label
    c     write (06,*) 'PLOTBOND nydd,nhbplot,nhbtot,nbresplot=',
    c    -  nydd,nhbplot,nhbtot,nbresplot
    c     print *,'PLOTBOND nclust,rdclust,numres,maxresdist=',
    c    -                  nclust,rdclust,numres,maxresdist
    c     print *,'PLOTBOND nframetot=',nframetot,' xtrajmax=',xtrajmax
    c      write (40,8711) (ihbtores(i),i=1,nhbtot)
    c8711  format(' IHBTORES:',/,(15i5))
            call roundlimint(nfrtot,ixd,nxd)
    c       Round up nframetot based on xtrajmax
            call roundlim(frtot,xdiv,nxd)
    c     xfac=xmm/float(nxd*ixd)
              call trnsfi(itemp1,index,nhbtot)
              call indexit(index,1,nhbtot,0)
          call roundlimint(ntotplot,iyd,nyd)
          call zeroiti(ihbon,0,nhbtot)
          call zeroiti(ihbstart,0,nhbtot)
            call readbitc(ires(1,ifr),ihb,nhbtot,30,MAXITEMS)
    c           ibx=index(ib)
              call trnsfi(ihbon,ihb,ntotplot)
    c         write (ips,1001) -xm*1.1,0.1*ym,' translate'
    c         Drawing bounding box
    c         write (ips,1017) trajnam(1:ltrajnam),ifirst,ilast,incr_traj
              call write_traj_lim(ips,' ',1,1,incr_tr,1)
              call psshow(ips,'Title: ',7)
              call psshow(ips,itit,ntit)
              call rgbcolor(ips,9)
    c         write (ips,1001) xm0+0.45*xmm,ym0-0.09*ymm09,' m'
              call psshow(ips,xlab,nxlab)
    c           Don't skip y label to leave room for res-res info
                call psshow(ips,ylab,nylab)
              call rgbcolor(ips,9)
    c             if (xdiv .ge. 100.0) then
    c           Print info
                  call makebondlab(iy,iy,iy-1,ibondlab,bondlab,lbondlab,
    c           Plot cluster marker bars
    c               if (mod(ncp,2) .eq. 0) xinc=1.01
    c         write (ips,*) 2,' setlinewidth'
    c       write (ips,*) lw,' setlinewidth'
                  call rgbcolor(ips,9)
    c           Draw line between ihbon(ibx) and ifr
                call rgbcolor(ips,-icolrb(ianc_anc(ibx)+1))
    c     Print tick on the right y axis
          call rgbcolor(ips,9)
    c     Print color code
            call rgbcolor(ips,-icolrb(1))
    c       write (ips,1001) xm0+10.0,ym0-0.10*ymm09,' m'
            call rgbcolor(ips,-icolrb(2))
          call rgbcolor(ips,9)
    c1017 format('(Trajectory analyzed:',a,' Frames',i7,' - ',i6,
    c    -  ' Increment:',i5,') show')
          character*8 atnames(maxrec),resnames(maxrsd)
          character*37 bondlab(maxbondlab)
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
          common /bondpairs/ ihbpair(2,MAXBONDS),ihb_pair_res(3,MAXBONDS)
    c     print *,'MAKEBONDLAB ilab1,ilab2,maxbondlab,ilab=',
    c    -  ilab1,ilab2,maxbondlab,ilab
            call blankout(bondlab(ilab-inc),1,37)
    c         res-aggregated plot
    c         write (6,9871) ilab,irr,ixr1,ixr2,ir1,ir2
    c9871     format(' ILAB,IRR=',2i4,' IXR1,IXR2=',2i4,' IR1,IR2=',2i4)
              call lastchar(resnames(ixr1),lc1,8)
              call lastchar(resnames(ixr2),lc2,8)
    c         Single bond plot
              call lastchar(resnames(ixr1),lc1,8)
              call lastchar(resnames(ixr2),lc2,8)
          character*(*) title,remark,xlab,y1lab,y2lab
    c     iyinc: increment the first index of y12 with iyinc (only for nplot=1)
    c     print *,'PLOT2FUN n=',n,' nplot=',nplot,' yl1=',y1lab(1:ly1lab)
    c     print *,'y2mn,y2dv,ny2dv=',y2mn,y2dv,ny2dv
    c     print *,'XLAB  =',xlab  (1:lxlab  )
    c     print *,'PLOT2FUN NOSD=',nosd
    c     do i=1,n
    c        write (40,9789) x(i),y12(1,i),sd12(1,i),y12(2,i),sd12(2,i)
    c9789    format(f6.0,' y1=',f9.4,' sd1=',f9.4,' y2=',f9.4,' sd2=',f9.4)
    c     end do
          call psheader(iplot,title,ntit,0,0,612,792,npspages,ipspage)
          call psshow(iplot,xlab,lxlab)
          call psshow(iplot,title,ntit)
            call write_traj_lim(iplot,' ',1,itrajfile,incr_tr,1)
    c     write (iplot,1002) ix0+ixwid,iy0+iyhgt+10
    c     Find X range
            call roundlim(xmax,xdiv,nx)
    c     Find Y1 range
          call arminmax2(y12,1,n,nplot,y1min,y1max,y2min,y2max,iyinc,
    c     print *,'y1min,y1max=',y1min,y1max,' ny1dv=',ny1dv
    c     print *,'y2min,y2max=',y2min,y2max,' ny2dv=',ny2dv
    c     print *,'IYINC=',iyinc
    c       Use 2nd column of y12 only
    c     Draw ticks, write axis values
    c       Draw ticks, write axis values
    c     Draw y1 tics
          call rgbcolor(iplot,-icol1)
    c     Print y1 axis tick values
          call rgbcolor(iplot,9)
            call rgbcolor(iplot,-icol2)
            call rgbcolor(iplot,9)
    c     if (xmin .ne. 0.0) then
    c       Plot initial value too
    c       ixmin=0
    c     else
    c       ixmin=1
    c     end if
          call rgbcolor(iplot,-icol1)
          call psshow(iplot,':',1)
    c       Print remark (if any)
            call rgbcolor(iplot,9)
            call psshow(iplot,remark,lremark)
            call rgbcolor(iplot,-icol1)
    c     Plot graphs
    c           Just draw a circle
    c           Plot error bars
            call rgbcolor(iplot,-icol2)
            call psshow(iplot,':',1)
            call psshow(iplot,y2lab,ly2lab)
    c         write (6,6734) i,x(i),xmin,xmax,y12(2,i),y2min,y2max
    c6734     format(i4,' X=',f8.3,' XMIN/MAX=',2f8.3,' Y12=',f8.3,
    c    -      ' Y2MIN/MAX=',2f8.3)
    c             Just draw a circle
    c             Plot error bars
            call rgbcolor(iplot,9)
    c1099  format('%%Trailer')
          character*(*) title
          call psshow(ips,title,lentit)
    c#    MMC routine 354 lstmod: 04/05/08
    c*****Postscript plot of n functions of the same variable
          character*(*) xlab,fclab(nf),tit,tit2
          character*1 marks(9)
    c     character*80 ident
    c     common /title/ ident(2)
    c     nf: number of functions to be plotted
    c     x: the x coordinates of the functions to be plotted
    c     ifg(if): First value in x for the if-th function
    c     imf(if), iml(if): Y {imf(if) - iml(if)} the value of the if-th function
    c     r0,cx: x coordonate labels are transformed as r0+x(i)*cx
    c     y00,yd: y scale minimum and unit, yd=0 => program finds them
    c     iprt: if .ne. 0, print the function values;
    c     tit: string containing the title; ntit: number of chars in tit
    c     print *,'ny,iprt,ntit,ntit2=',ny,iprt,ntit,ntit2
    c     print *,'nxmax,nymax,nf=',nxmax,nymax,nf
    c     print *,'tit=',tit(1:ntit)
    c     print *,'tit2=',tit2(1:ntit2)
    c     print *,'fclabs=',(fclab(i)(1:lfclab(i)),i=1,nf)
    c     print *,'PLOTNPS ipspage,npspages,iout=',ipspage,npspages,iout
          call psheader(iplot,tit,ntit,0,0,612,792,npspages,ipspage)
    c     write (iplot,1002) ix0+10,iytit
    c     write (iplot,1004) ident(1)
    c     iytit=iytit-15
    c     write (iplot,1002) ix0+10,iytit
    c     write (iplot,1004) ident(2)
    c     write (iplot,1010) ix0,iy0,ixwid,iyhgt,-ixwid
    c     Drow graph boundary box
    c     Draw ticks, write axis values
    c       Mark residues
    c       write (6,9292) (imarks(i),i=1,120 )
    c9292    format(50i1)
    c     Plot graphs
    c     iyhgt=480
    c       print *,'if,imf(if),iml(if)=',if,imf(if),iml(if)
    c         write (77,*) 'if,ig,y(ig)=',if,ig,y(ig)
    c1010  format('% Drawing of graph boundaries',/,'newpath',/,
    c     -  i3,1x,i3,' moveto',/,i4,' 000 rlineto',/,'000 ',i4,' rlineto',/,
    c     -  i5,' 000 rlineto',/,'closepath',/,'stroke')
    c#    MMC routine 355 lstmod: 12/05/03
    c*****Set the PS color to frac way on the rainbow scale (0<=frac<=1)
          character*(*) title,plotdesc,remark,xlab,ylab,timelab
    c     print *,'PLOT2D xmn,xdv,nxdv=',xmn,xdv,nxdv
    c     print *,'PLOT2D ymn,ydv,nydv=',ymn,ydv,nydv
    c     print *,'START PLOT2D iout,npspages,ipspage=',
    c    - iout,npspages,ipspage
          call psheader(iplot,title,ntit,0,0,612,792,npspages,ipspage)
          call psshow(iplot,xlab,lxlab)
          call psshow(iplot,ylab,lylab)
          call psshow(iplot,title,ntit)
    c     write (iplot,1002) ix0+ixwid,iy0+iyhgt+10
    c       Print remark (if any)
            call psshow(iplot,plotdesc,lplotdesc)
    c       Print remark (if any)
            call psshow(iplot,remark,lremark)
    c     Find X,Y range
          call arminmax2(xy,1,n,2,xmin,xmax,ymin,ymax,0,2)
    c     Draw ticks, write axis values
    c       Plot initial value too
    c     Draw ticks, write axis values
    c       Plot initial value too
    c       Left tick
    c       write (iplot,1011)
    c       Right tick
    c       write (iplot,1012)
    c       Plot initial value too
    c       Lower tick
    c       write (iplot,1011)
    c       Upper tick
    c       write (iplot,1012)
    c     Plot traces
          call rrgbcolor(iplot,1,n,1)
    c             Don't draw lines for moves under periodicity
    c               Draw 2-color line
                    call rrgbcolor(iplot,i,n,1)
                  call rrgbcolor(iplot,i,n,1)
          call rainbowscale(iplot,ix0,ixwid,iy0-100,n,0.0,0.0,0.0,
    c1099  format('%%Trailer')
    c     if iyinc=1, use the 2nd column of ar (valid only for nplot=1)
    c     print *,'ARMINMAX2 nfst,nlst=',nfst,nlst
    c     For dmaxinp, find rounded div and ndiv such that div*ndiv ~ dmaxinp
    c     For nmaxinp, find rounded idiv and ndiv such that idiv*ndiv ~ nmaxinp
    c     print *,'ROUND idiv,ifac,nmaxinp=',idiv,ifac,nmaxinp
    c     print *,'ROUND2 idiv,ndiv=',idiv,ndiv
    c      print *,'armin(inp)=',armin
    c       First, round off the minimum
            call roundlim(armin,xymindiv,nxymindiv)
          call roundlim(range,xydiv,nxydiv)
    c     print *,'CHECKBEND n=',n
          call fitpoints(c,n,3,x0,rn,axfact,idebughx)
              call writepdbd(77,c(1,i),i,i,'C   ','CA ','A',1.0,0.0)
            call writepdbd(77,x0,n+1,n+1,'S   ','X0 ','B',1.0,0.0)
              call dvsum(x0,rn,x1)
              call writepdbd(77,x1,n+2,n+2,'S   ','RN ','B',1.0,0.0)
    c     Find the shortest axis-Calpha distance
    c     "Pull" uniformly all the Calphas toward the axis as close as possible
              camod(k,i)=c(k,i)+rr*perpvec(k,i)
              call writepdbd(77,camod(1,i),i,i,'C   ','CA ','B',1.0,0.0)
    c     Calculate the radius of the fitting circle to the pulled points
          call circfit(c,n,axisdir,circ)
    c     Find the projection of the pulled points on the fitting plane
            call dvdif(camod(1,i),x0,x1)
              camod(k,i)=camod(k,i)-rr*rn(k)
          call fitpoints(camod,n,2,x0,rm,axfact,0)
              call writepdbd(77,camod(1,i),i,i,'N   ','CAP','C',1.0,0.0)
            call writepdbd(77,x0,n+1,n+1,'P   ','X0 ','B',1.0,0.0)
            call dvsum(x0,rm,x1)
            call writepdbd(77,x1,n+2,n+2,'P   ','RM ','B',1.0,0.0)
    c       Check if projections are indeed in a plane
            call dvdif(camod(1,1),x0,x1)
            call dvdif(camod(1,2),x0,x2)
            call dvprd(x1,x2,x)
              call dvdif(camod(1,i),x0,x1)
            call trnsfrd(cx,camod,3*n)
    c     Get the distance vectors between the point and its projection on the line
              camod(k,i)=(camod(k,i)-x0(k))-axfact(i)*rm(k)
              call dvdif(cx(1,i),camod(1,i),x)
              call writepdbd(77,x,i,i,'O   ','CA0','D',1.0,0.0)
    c       Eliminate (and count) residuea that are within tolerance of the axis
                call trnsfrd(camod(1,i-nnear),camod(1,i),3)
          call zeroiti(nupdown,0,2)
          call runtest(nup,ndown,nrun,idecide)
    c*****Computes test for correlation, type of correlation
    c     Minimum critical values for correlation test:
    c     Maximum critical values for correlation test:
    c       High
    c       Low
    c       Low
    c       Correlated
    c       Correlated
    c       Uncorrelated
    c     For n points in c, find the mean (c0) and the direction of the normal
    c     to the best fitting plane (ndim=3) or the direction of the best fitting
    c     line (ndim=2)
    c     print *,'FITPOINTS n=',n
          call zeroitd(c0,3)
              c0(k)=c0(k)+c(k,i)
            c0(k)=c0(k)/dfloat(n)
    c     Find the eigenvectors a and eigenvalues mu of rr
          call dtred2(rr,3,3,diag,offdiag)
          call dtqli(diag,offdiag,3,3,rr,ierr)
    c     The columns of the matrix rr are the eigenvectors
    c     The eigenvalues are in diag
            call dvdif(c(1,i),c0,dd)
          character*(*) label
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
    c     print *,'FINDBESTREP iw0,icl,nframe=',iw0,icl,nframe
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
          character*80 title
          character*(*) title2
          character*(*) xlab
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
          character*4 yclab(1)
          character*80 title2add
    c     print *,'PLOT2DRMSD isymm,nframe,nframeref=',
    c    -                   isymm,nframe,nframeref
    c     print *,'PLOT2DRMSD noclose,matplotonly=',noclose,matplotonly
    c     print *,'PLOT2DRMSD ioclose,matplotonly=',noclose,matplotonly
          cdiammin=100000.0
          cdiammax=0.0
            call indexit(indexa,1,nframe,0)
            call findbestrep(iw0,0,0,nframe,indexa,irdavmn,irdmxmn,
    c       print *,'IYDEL,IYTOP,IYRANGE,NFRAMEREF=',
    c    -    iydel,iytop,iyrange,nframeref
            call contractmat(rmsd2d,nframe,nframeref,nframeplot,
    c       iymax=iydel+iyrange
    c       iymax=iymax+15
              call psshow(iw1,
            call psshow(iw1,title2add(1:ladd),ladd)
              call indexit(ixshuffle,1,nframe,0)
              call indexit(ixshuffleref,1,nframeref,0)
            call plotmat(iw1,kc,rmsd2d,dc,nframeplot,nframerefplot,0,0,0,0,
              close (iouttemp,status='delete')
            call colcodeminmax(iw1,20+ixcent,-iydel,nrep,ncolcode,
    c       Draw lines in the matrix delineating the clusters
    c          print *,'NCLX,NCLY=',nclx,ncly,' NFRAME=',nframe
    c          write (6,9876) (ilastclx(i),i=1,nclx)
    c9876      format('ILASTCLX:',/,(20i4))
              call rgbcolor(iw1,9)
    c         ixdell=ixdelsh
    c    -      ixdell=ixdelsh+40/scalefac
    c       2D RMSD map, dont calculate normalization
              call arminmax2(res(1,1,8),1,nframe-1,2,armin1,armax1,armin2,
              call roundlim(armax1,y1div,ny1div)
              call roundlim(armax2,y2div,ny2div)
              call plot2fun(iw1,2,xtraj,res(1,1,8),res(1,1,9),nframe-1,
          call zeroiti(nrmsd,0,100)
    c       dcoef=rcmax/sqrt(float(nframe))
              call roundlim(pmax1,y1div,ny1div)
              call roundlim(pmax2,y2div,ny2div)
              call plot2fun(iw1,2,xrmsd,res(1,1,9),res(1,1,9),100,0.0,
    c       Cross RMSD map
            call roundlim(pmax1,y1div,ny1div)
            call plot2fun(iw1,1,xrmsd,res(1,1,9),res(1,1,9),100,0.0,
          character*132 line(maxrec)
          character*80 title
          character*(*) trajfile
          common /nnwork/ trajdist(MAX2D,MAX2D),cav(3,MAX2D),
          common /colorinfo/ ncolcode,maxcolcode
          character*4 yclab(1)
          character*17 normtyp(4)
    c     print *,'PLOT ATOMDIST_SD nframe,ndist=',nframe,ndist
              cav(k,ia)=cav(k,ia)/dfloat(nframe)
              cavsng(k,ia)=cav(k,ia)
            cavs(ia)=dsqrt(cavs(ia)/dfloat(nframe)-rr)
          call plotmat(iw1,kc,rc,trajdist,ndist,ndist,0,0,0,0,navg,
          call colcodeminmax(iw1,20+ixcent,-iydel,nrep,ncolcode,
          character*500 label(mx2d)
          character*500 line
          character*(*) inpfile
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
          call openfile(iw0,0,'distance matrix file',20,'old',inpfile,
          call askyn('Do you want to limit the label to its last word',47,
          call checkdim(n,MAX2D,'MAX2D',5,'matrix dimension',16,0)
            call blankout(line,1,llen)
            call lastchar(line,lline,llen)
            call nextblank(line,ic,llen)
    c       print *,'ic1,ic=',ic1,ic
    c       print *,'line=',line(1:30)
            call blankout(label(i),1,llen)
              call nextchar(line,ic,llen)
              call laststring(line,icf,icl,lline,llen)
    c       print *,'i=',i,' ix=',ix(i),' label=',label(i)(1:lline-ic+1)
          close (iw0)
    c     do i=1,n
    c       print *, ix(i),label(i)
    c       print *, (rmsd2d(i,j),j=1,n)
    c     end do
          call transform_mat(rmsd2d,n,mx2d,rijrmin,iuout)
          call transform_dist(rijrmin,rijmin)
          call transform_dist(rijrmax,rijmax)
          call zeroiti(n10,0,10)
            call askyn('Do you want a list of identical pairs',37,
                    call lastchar(label(i),lci,500)
                    call lastchar(label(j),lcj,500)
          character*1 ans
          common /transform_dist_dat/ itranstyp,nexp,rijmin,fracexp
          call quiz(ans,itranstyp,'u',' ',0,'transformation type',19,
              call getint('Exponent',8,1,1,6,nexp,0)
                call transform_dist(r(i,j),r(i,j))
          common /transform_dist_dat/ itranstyp,nexp,rijmin,fracexp
          character*(*) trajname,trajname2,system
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
          character*16 errlab
          character*80 line
          call blankout(errlab,1,16)
            call blankout(line,1,80)
    c         RMSD data found
    cConfig #     1 (inp #:          1) Ref #     1 (inp #:         1)
    c         read (line(52:57),*,err=888) jx
                call blankout(line,1,80)
                call lastchar(line,lc,80)
              call lastchar(line,icl,80)
              call lastchar(line,icl,80)
              call save_traj_lim(ifirst,ilast,incr,1)
              call save_traj_lim(ifirst2,ilast2,incr2,2)
              call lastchar(line,icl,80)
              call blankout(line,1,80)
              call lastchar(line,lsystem,80)
    c       Cross RMSD map
    c     write (78,2001) nframex,nframey
    c     do i=1,10
    c       write (78,*) 'RMSFCLUSTER i=',i,'RMSD:'
    c       write (78,8712) (rmsd2d(i,j),j=1,nframex)
    c8712   format(10f8.4)
    c     end do
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
          character*6 lab
          call zeroiti(nsame1,0,nframe1)
          call zeroiti(nsame2,0,nframe2)
    c     print *,'MATCHTRAJ iw0,nframe,nframeref=',iw0,nframe,nframeref
    c       Find closest to frame if
    c       Find closest to frame if
    c     Calculate the RMSD of the two structures after obtaining the best fit
    c     using Kabsch formula
          character*2 iatnm2
          common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
    c     print *,'Start RMSD nframeref,nframe=',nframeref,nframe,
    c    -  ' nfinalov,nfinalrmsd=',nfinalov,nfinalrmsd
    c     Calculate rMSD w/o overlay
    c      write (iout,8811) n,nfinalrmsd,nfinalov,iedit,limresrange
    c8811  format(' n,nfinalrmsd,nfinalov,iedit,limresrange=',3i8,2i3)
    c      write (iout,8812) 'indxov1',(indxov1(i),i=1,nfinalov)
    c      write (iout,8812) 'indxrmsd1',(indxov2(i),i=1,nfinalrmsd)
    c8812  format(1x,a,':',/,(20i4))
    c      write (iout,8813) 'C',((c(k,indxov1(i)),k=1,3),i=1,nfinalov)
    c      write (iout,8813) 'CO',((co(k,indxov2(i)),k=1,3),i=1,nfinalov)
    c8813  format(1x,a,':',/,(5f10.5))
    c      write (iout,8814) (atw(indxov1(i)),i=1,nfinalrmsd)
    c8814  format(' ATW:',/,(5f10.5))
            call bestoverlay(nfinalov,indxov1,indxov2,co,c,atw,atwsum,c1,c2,
    c       Shift the full set of coordinate to com1,com2, and rotate by rot
            call shiftmol(co,n,com1,c1,-1.0)
            call shiftmol(c,n,com2,c2,-1.0)
            call rotate_c(c2,n,rot,c2,'RMSD',4)
    c      do i=1,45
    c        write (78,8674)i,(co(k,i),k=1,3),(c1(k,i),k=1,3),(c2(k,i),k=1,3)
    c8674    format(i3,' co=',3f10.5,' c1=',3f10.5,' c2=',3f10.5)
    c      end do
          call trajlimtest(nframe,MAXFRAMES)
    c     print *,'RMSF nslt=',nslt
    c      write (78,9878) (i,index(i),index(i),atw(index(i)),
    c     -  (cref(k,index(i)),k=1,3),(c(k,index(i)),k=1,3),i=1,nslt)
    c9878  format(i5,' ix1,2=',2i4,' aw=',f8.3,' c1=',3f10.5,' c2=',3f10.5)
    c     write (40,8789) (index(ia),ia=1,nslt)
    c8789 format(20i4)
    c         write (76,9782) ia,iaa,ir,atw(iaa),(cref(k,iaa),k=1,3)
    c9782     format(' ia,index(ia),ir=',3i5,' aw=',f8.3,' cref=',3f10.5)
    c       bfacavg(ir)=bfacavg(ir)+sqrt(rmsfs/atwrsum)
          call zeroitd(rmsfsum,nresslt)
    c          write (77,9671) ia,ir,index(ia),atw(index(ia)),
    c     -      (cdp(k,ia),k=1,3),(cdp2(k,ia),k=1,3)
    c9671      format(i4,' ir=',i3,' ix=',i4,' atw=',f8.3,' cdp=',3e12.5,
    c     -      ' cdp2=',3e12.5)
    c       write (77,*) 'ir=',ir,' RMSFAV=',rmsfav(ir)
    c     Calculate the COM of the solute
    c     Calculate the solute and the whole cell dipole moment when icharges >0
          common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
          call trajlimtest(nframe,MAXFRAMES)
          call zeroitd(dipole,3)
          call zeroitd(com,3)
                com(k)=com(k)+c(k,index(ia))*aw(index(ia))
            com(k)=com(k)/awsum
          call zeroitd(dipole,3)
          character*(*) label
          character*(*) label2d(mx2d)
          character*18 memlab
          character*200 memdir,memfilename
          character*500 line
          character*1 ans
          character*41 clstyp
          common /cluster_typ/ nclstyp,inumclst(9),ireadcutoff(9),
    c     Clustering
    c     print *,'CLUSTERD ndim,rdmin,rdmax,ietotsaved=',
    c    -  ndim,rdmin,rdmax,ietotsaved
    c      write (iout,7934) (index2d(i),i=1,ndim)
    c7934  format(' CLUSTERDISTR INDEX2D=',15i4)
              call askyn('No clustering was done - still want to quit',43,
            call getint('Number of clusters requested',28,999999,1,ndim,
            call getreal(line,22+llabel,rdclustdef,rdclust,1,ihelp)
    c       Set density type and neighbor min
            call quiz(ans,idenclstyp,' ',' ',0,'density clustering variant',
          call zeroit(rdlim,ndim)
          call trnsfi(it2d,index2d,ndim)
    c       Just cluster using rdclust as the threshold
            call rmsdcluster(rdclust,1,ndim,index2d,iwt,ixclst,ifclst,
    c      write (iout,7824) (ifclst(i),ilclst(i),i=1,nrdclust)
    c7824  format(' After RMSDCLUSTER Cluster limits: ',('[',i5,',',i5,']'))
    c      write (iout,6792) 'IXCLST',(ixclst(i),i=1,ndim)
    c6792  format(' AFTER RMSDCLUSTER ',a,':',/,(20i5))
    c       Members of cluster ic: (index2d(i),i=ifclst(ic),ilclst(ic))
            call reportclust(ndim,0,1,nclust,ifclst,ilclst,index2d,value,
    c       Vary the threshold until nrdclust clusters result
            call zeroiti(indexa,0,ndim)
            call trnsfi(it4,index2d,ndim)
              call trnsfi(index2d,it4,ndim)
              call rmsdcluster(rdclust,1,ndim,index2d,iwt,ixclst,ifclst,
            call reportclust(ndim,0,1,nclust,ifclst,ilclst,index2d,value,
            call askyn('Do you want uniformly spaced cutoffs',36,1,1,
              call transform_dist(rdmin,rdtmin)
              call transform_dist(rdmax,rdtmax)
                call getreal('Largest cutoff',14,rdtmax,rcutulim,1,00)
                call getreal('Smallest cutoff',15,rdtmin,rcutllim,1,00)
                call getreal('Largest cutoff',14,rdtmin,rcutulim,1,00)
                call getreal('Smallest cutoff',15,rdtmax,rcutllim,1,00)
                cutofflist(icut)=
                call getreal(line,11,999999.0,cutofflist(icut),1,00)
            call getint('Cluster level to calculate average',34,ncutoff,1,
            call askyn('Do you want to write cluster member files',41,1,-1,
              call checkdir(memdir,lmemdir,iu_clst,iopen)
            call indexit(index2d,1,ndim,0)
    c           write (iout,*)
    c    -        'cutofflist(icut),ifcl_prev(icl),ilcl_prev(icl)=',
    c    -        cutofflist(icut),ifcl_prev(icl),ilcl_prev(icl)
                call rmsdcluster(cutofflist(icut),ifcl_prev(icl),
              call trnsfi(ifcl_prev,ifclst,ncltot)
              call trnsfi(ilcl_prev,ilclst,ncltot)
    c           Find out if cluster ends at level lev_avg
                call lastchar(label2d(index2d(ia)),lc,llabel2d)
    c           Write member list file
                  call writeint(memfilename,lmemfilename+1,
    c           call openfile(iu_clst,0,' ',1,'new',memfilename,
    c    -        lmemfilename,notfound,0,1,1,1,0)
    c             Write list
                    call lastchar(label2d(index2d(imem)),lc,llabel2d)
                  close (iu_clst)
                call laststring(label2d(index2d(ia_rep)),ifc,ilc,lc,500)
    c       Sub-clustering only works for single-link clustering
            call askyn('Do you want to try sub clustering',33,1,-1,isubcl,
              call quiz(ans,isubclustertyp,' ',' ',0,
              call getint('Minimum number of members for subclustering',43,
                call getreal(
                call getreal(
                  call rmsdsubcluster(ic,ifclst(ic),ilclst(ic),index2d,
    c               Adjust cluster limits
                    call sortlist(iout,index2d,ilclst(nclust),it1,it2,'IX2',
                    call reportclust(ndim,ic,ic,ic+nclustic-1,ifclst,ilclst,
                call lastchar(label2d(ixclst(i)),llab,llabel2d)
            call askyn(
            call askyn(
                  cv(i-ifclst(ic)+1)=etotsaved(1,index2d(i))
          call askyn('Do you want to run more clustering',34,1,-1,
            call trnsfi(index2d,it2d,ndim)
          character*(*) xtrajlab
    c     Plot the history of cluster membership
    c     print *,'CLUSTERPLOT ncl=',ncl,' nframe=',nframe,' ips=',ips
          call plot2fun(ips,1,xtraj,value,value,nframe,0.0,0.0,0,0.0,
          character*(*) label
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
    c     Set up neighbor list based on the threshold
    c     print *,'RMSDCLUSTER iout,n,iclstyp,MX2D=',iout,n,iclstyp,MAX2D
    c     print *,'RMSDCLUSTER iclstyp,idenclstyp,nnmin=',iclstyp,idenclstyp,nnmin
    c     write (40,*)'RMSDCLUSTER  nfrst,n,iclstyp=',nfrst,n,iclstyp
    c     write (40,*)'RMSDCLUSTER  RMSDCLUST=',rmsdclust
    c      write (40,7934) (index2d(i),i=1,n)
    c7934  format(' RMSDCLUSTER INDEX2D=',15i4)
    c        Cluster has one element
          call zeroiti(nng,0,n)
            call trnsfi(index2dd,index2d,n)
            call indexit(index2dd,1,n,0)
    c     write (77,*) 'RMSDCLUSTER nfrst,n=',nfrst,n
    c     write (79,*) 'nfrst,n=',nfrst,n
    c     do i=1,n
    c       write (79,1004) (rmsd2d(i,j),j=1,n)
    c     end do
    c1004 format(5e13.6)
    c         if (rmsd2d(i,j) .eq. 0.0) then
    c           write (6,*) 'ii,jj=',ii,jj,'i,j=',i,j
    c         end if
    c       Clique-clustering - sort lists
    c       Sort iorder by iwt
            call indexit(iorder,1,n,0)
            call mrgsrt(6,iorder(nfrst),t2,nmem,it2,it3,it4,t1,nmem)
    c       Sort neighbors by iwt
                call mrgsrt(6,ing(1,i),t2,nng(i),it2,it3,it4,t1,nng(i))
    c      do ii=nfrst,n
    c        i=index2dd(ii)
    c        write (40,8967) ii,i,(rmsd2d(i,index2dd(j)),j=1,n)
    c8967    format(' ii,i=',2i5,' RMSD2D:',15f5.2,/,(15f5.2))
    c        write (40,8968) ii,i,(ing(jj,ii),jj=1,nng(ii))
    c8968    format(' ii,i=',2i5,' ING:',15i4,/,(14i4))
    c      end do
    c      do i=1,n
    c        write (6,8701) i,nng(i),(ing(j,i),j=1,nng(i))
    c8701    format(i4,' nn=',i4,' in=',10i5)
    c      end do
    c       Single-link clustering
            call clstrs(ing,nng,it1,nfrst,n,iclst,ifirst,ilast,0,nofcls,
    c       K-medoids clustering - no COM
            call clstrs_kmedoids(nfrst,n,index2dd,iclst,ifirst,ilast,nofcls,
    c         Extract coordinates of center nodes
                call trnsfr(cent(1,ic),c(1,icent_fin(ic)),3)
    c       K-means clustering -  COM-based
            call clstrs_kmeans(nfrst,n,index2dd,iclst,ifirst,ilast,nofcls,
            call clstrs_maxnn(ing,nng,nfrst,n,iclst,ifirst,ilast,0,
    c       Clique-based clustering
    c       iout_c=6
            call clstrs(ing,nng,it1,nfrst,n,iclst,ifirst,ilast,0,nofcls,
            call clstrs_density(nfrst,n,iclst,ifirst,ilast,nofcls,nofcls_t,
    c     Sort iclslt within each cluster
    c      write (06,9671) 'U',nofcls,(iclst(i),i=1,n)
    c9671  format(1x,a,' NCLUST=',i4,' ICLST:',/,(20i4))
              call indexit(it1,1,nmem,0)
              call mrgsrt(6,it1,t2,nmem,it2,it3,it4,t1,nmem)
    c     Rearrange the elements of index2d in the order of iclst
          call trnsfi(index2d(nfrst),it1(nfrst),n-nfrst+1)
          character*(*) label
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
    c     Break up a cluster by removing low-neighbor number nodes
    c     First recreate the nn list for the cluster selected in a way
    c     to stay within the cluster
          call zeroiti(nng,0,n)
          call clstrs(ing,nng,ix1,nfrst,n,iclst,ifirst,ilast,0,nofcls,i2,0,
          call trnsfi(i3,nng,n)
    c         Subclustering by # of neighbors
              call trnsfi(i4,nng,n)
    c         Find NN range (singletons excluded)
    c         Now eliminate all nodes with nnmin neighbors
    c         Subclustering by 'density'
    c       i5 has the list of nodes to delete
            call clstrs(ing,nng,i1,nfrst,n,iclst,ifirst,ilast,0,nofcls,
    c       Count the number of non-singleton clusters
            call askyn('Do you want to shave more nodes',31,
    c       Add back to the real clusters the deleted nodes
            call zeroiti(i5,0,n)
    c       Now i5 is the cluster number of a remaining node or zero
    c           Find a remaining neighbor, assign ia to its cluster
    c       Now, recluster the nodes based on i5
            call indexit(i1,1,n,0)
            call mrgsrt(6,i1,value,nmem,i2,i3,i5,t1,nmem)
    c       Finally, rearrange the elements of index2d in the order of iclst
            call trnsfi(index2d(nfrst),i1(nfrst),n-nfrst+1)
          character*(*) label
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
    c     print *,'REPORTCLUST nfclst,nlclst,isorttype,ifindbestrep=',
    c    -  nfclst,nlclst,isorttype,ifindbestrep
    c      write (77,8822) (index2d(i),i=1,ilclst(nlclst))
    c8822  format(' index2d:',25i5)
    c     do i=1,ndim
    c       write (77,9782) i,(rmsd2d(i,j),j=1,ndim)
    c9782   format(i4,' ccc=',30f5.2)
    c     end do
          call indexit(ixrank,1,nclust,0)
          call mrgsrt(6,ixrank,rankav(nfclst),nclust,ifa_s,ila_s,ih,cv,
          cdiammin=100000.0
          cdiammax=0.0
            call indexit(indexa,1,ndim,0)
            call findbestrep(0,0,0,ndim,indexa,irdavmn,irdmxmn,
    c         Sort by increasing member index
              call indexit(it1,1,nmem,0)
              call mrgsrt(6,it1,value,nmem,ifa_s,ila_s,ih,cv,max2d)
    c         Sort by decreasing occurrence
              call mrgsrt(6,it1,value,nmem,ifa_s,ila_s,ih,cv,max2d)
              call trnsfi(index2d(ifclst(ic)),it1,nmem)
              call zeroiti(indexa,0,ndim)
                call findbestrep(iout,ic,0,ndim,indexa,irepav(ic),
                call findbestrep(iout,icl0,ic,ndim,indexa,irepav(ic),
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
    c     Count the  number of pairs in each cluster that are within rmsdsim
          call getreal('MAXimum RMSD for intra-cluster similarity',41,
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
    c     Count the  number of pairs in each cluster that are within rmsdsim
          character*(*) trajnam1,trajnam2
          common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
    c     Map traj 2 ont clusters of traj1
    c     print *,'MAPCLUSTX nclust1,nclust2,iout=',nclust1,nclust2,iout
          call zeroiti(nndistsum,0,nclust1)
          common /nnwork/ trajcorr(MAX2D,MAX2D),cav1(3,MAX2D),
    c     print *,'RESIDCORRR nframe,ncorr=',nframe,ncorr
            call zeroitd(cav1,3*ncorr)
            call zeroitd(cav2,3*ncorr)
            call zeroitd(cavs1,ncorr)
            call zeroitd(cavs2,ncorr)
              call zeroitd(trajcorr(1,ir),ncorr)
              cav1(k,ir)=cav1(k,ir)+c1(k,index(ir))
              cav2(k,ir)=cav2(k,ir)+c2(k,index(ir))
              cavs1(ir)=cavs1(ir)+scprod(c1(k,index(ir)),c1(k,index(ir)))
              cavs2(ir)=cavs2(ir)+scprod(c2(k,index(ir)),c2(k,index(ir)))
          common /nnwork/ trajdist(MAX2D,MAX2D),cav(3,MAX2D),
    c     print *,'ATOMDIST_SD nframe,ndist=',nframe,ndist
            call zeroitd(cav,3*ndist)
            call zeroitd(cavs,ndist)
              call zeroitd(trajdist(1,ir),ndist)
              cav(k,ir)=cav(k,ir)+c(k,ldist(ir))
              cavs(ir)=cavs(ir)+scprod(c(1,ldist(ir)),c(1,ldist(ir)))
          common /nnwork/ ccc(MAX2D,MAX2D),itemp(MAX2D),fill(IFILL2)
    c     Calculate the normal modes from the residue correlation matrix
          common /nnwork/ trajcorr(MAX2D,MAX2D),drow(MAX2D),
          character*80 lineinp
    c     print *,'NORMALMODES inpt=',inpt,' inptyp=',inptyp
    c     Read in the covariance matrix
    c       Read matrix from the temp file written by Simulaid
    c       Binary
    c       Ascii
            call lastchar(lineinp,lc,80)
    c     write (78,*) 'ncorr=',ncorr
    c     do i=1,ncorr
    c       write (78,1004) (trajcorr(i,j),j=1,ncorr)
    c     end do
          call dtred2(trajcorr,ncorr,MAX2D,diag,offdiag)
          call dtqli(diag,offdiag,ncorr,MAX2D,trajcorr,ierr)
    c     Print normal modes
    c     Columns of trajcorr are the eigenvectors
          call askyn('Do you want to sort by eigenvalues',34,1,-1,isortev,0,
          call indexit(index,1,ncorr,0)
            call mrgsrt(6,index,value,ncorr,ifa,ila,itemp,temp,maxt)
              call write_traj_lim(iout,
          character*(*) title
          common /nnwork/ trajcorr(MAX2D,MAX2D),cav1(3,MAX2D),
          character*4 yclab(1)
          call indexit(ixshuffle,1,MAX2D,0)
          call write_traj_lim(iout,'Residue covariances and correlations',
              cav1(k,ir)=cav1(k,ir)/nframe
              cav2(k,ir)=cav2(k,ir)/nframe
            cavs1(ir)=cavs1(ir)/nframe
            cavs2(ir)=cavs2(ir)/nframe
    c     Save the covariance matrix for eigenvalue calculation
            call openfile(iucorrmat,0,'log',3,'new','trajcorr.mat',12,
            call plotmat(iplt,mx,rx,trajcorr,ncorr,ncorr,0,0,0,0,1,nrep,30,
            call rainbowscale(iplt,ixdelsh+50,450,iydel,0,0.0,rcmin,rcmax,
            call plothead(iplt,xm,ym-15,title,0,
            cavs1(ir)=dsqrt(dabs(trajcorr(ir,ir)))
          call plotmat(iplt,mx,rx,trajcorr,ncorr,ncorr,0,0,0,0,1,nrep,30,
          call rainbowscale(iplt,ixdelsh+50,450,iydel,0,0.0,-1.0,1.0,
          call plothead(iplt,xm,ym-15,title,0,
    c#    MMC routine 462 lstmod: 07/15/02; clicque option added
    c*****Find all clusters in a network and sort atoms in a cluster by groups
    c     Input parameters:
    c     n0,n: Use atomnumbers (vertices) from n0 to n (inclusive)
    c     maxat,maxneig,maxgr: Array sizes - see dimension statement above
    c     nneig(i) : number of neighbours (bonded) of molecule i
    c     ineig(j,i) : j-th neighbour of atom i
    c     Workspace arrays:
    c     iused(i) : 1 - atom i is not accounted for yet
    c                0 - atom i is already accounted for
    c     nnloop (i) : copy of nneig in loops 1 and 2 (temporary storage)
    c                  the number of loop-closing bonds of atom i thereafter
    c     Output parameters:
    c     nofcls0: Number of disconnected clusters (groups) found previously
    c     nofcls: Number of disconnected clusters (groups) found
    c     iclst,ifirst,ilast: The elements of the ig-th group are atoms
    c     (iclst(ia),ia=ifrst(ig),ilast(ig))
    c
    c     Description of the algorithm:
    c     Starting with an atom, the algorithm successively includes its
    c     neighbours and then the neighbours of the atoms already on the list.
    c     By excluding atoms that 're-occurred' the algorithm essentially generates
    c     a spanning tree.
    c     Initialization
    c     write (iout,*) 'CLSTRS n0,n,nofcls0,iverb=',n0,n,nofcls0,iverb
    c         Start search
    c         ncl is the number of elements in the cluster
    c         ic is the index of the atom under consideration
    c         kroot is the serial no of the lowest element in the cluster
    c         that may still have neighbours not examined yet
    c         NTOTSKIP=0
    c                 Include ic into the list
    c                 Check if ic is connected to all members found so far
    c                   No match - skip
    c                   New clique member found
    c                   NTOTSKIP=NTOTSKIP+1
    c                write (iout,8734) (iclst(ii),ii=nfound+1,nfound+ncl)
    c8734            format(' Clique so far:',(20i3))
    c               Now search for the first unused neighbor of ic
    c             if (NTOTSKIP .gt. 100) stop
    c           Neighbour chain ended, back to the kroot-th atom in the list
    c         Cluster of ncl elements found
    c         stop
    c         memmax=0
    c         memmin=10000000
    c         do ia=nfound+1,nfound+ncl
    c           if (memmax .lt. iclst(ia)) memmax=iclst(ia)
    c           if (memmin .gt. iclst(ia)) memmin=iclst(ia)
    c         end do
    c     ilast(nofcls)=n
    c     il=ifirst(nofcls)+ncl-1
    c*****Find all clusters in a network and sort atoms in a cluster by groups
          character*1 ans
    c     Input parameters:
    c     n0,n: Use atomnumbers (vertices) from n0 to n (inclusive)
    c     maxnode,maxgr: Array sizes - see dimension statement above
    c     Workspace arrays: it1,it2,it3
    c     Output parameters:
    c     nofcls: Number of clusters requested
    c     iclst,ifirst,ilast: The elements of the ig-th group are atoms
    c     (iclst(ia),ia=ifrst(ig),ilast(ig))
    c
    c     Description of the algorithm:
    c     k-means -> k-medoids
    c     Cluster centers have to be one of the nodes.
    c     The center node has the smallest max dev from the rest of the cluster
    c     print *,'CLSTRS_K start n0,n,nofcls,mx2d=',n0,n,nofcls,mx2d
    c     Initialization
    c      write (6,9671) (index2d(i),i=n0,n)
    c9671  format(' CLSTRS_KMEANS INDEX2D:',/,(20i4))
          call zeroiti(icent,0,nitersave*nofcls)
          call init_kmedoids(nofcls,index2d,ifirst,ilast,icent,rmsd2d,n0,n,
          call quiz(ans,icenttyp,'r',' ',0,'k-medoids center choice',23,0,
          call sortlist(0,icent,nofcls,it1,it2,'IK0',0,mx2d)
    c     print *,'CENT=',(icent(i),i=1,nofcls)
          call zeroiti(indxclst,n0-1,n)
    c       Find the nearest center to each node
    c       write (6,7711) 'Bef sort indxclst:',(indxclst(i),i=n0,n)
    c7711   format(' KMEANS ',a,(/i3,19i4))
    c         Find the new centers
              call indexit(it1,1,nnode,0)
    c         do i=n0,n
    c           w(i-n0+1)=indxclst(i)
    c         end do
              call mrgsrti(6,it1,indxclst(n0),nnode,it2,it3,it4,it5,n)
    c         write (6,7711) 'aft sort indxclst:',(it1(i),i=1,n-n0+1)
    c         write (6,7711) 'aft sort it1:',(it1(i),i=1,n-n0+1)
                    call mrgsrt(6,indxmax,distmax,nmem,it2,it3,it4,t1,nmem)
                    call mrgsrt(6,indxav,distav,nmem,it2,it3,it4,t1,nmem)
    c                write (77,9873) (indxmax(i),i=1,nmem)
    c                write (77,9874) (distmax(i),i=1,nmem)
    c                write (77,9875) (it2(i),i=1,nmem)
    c                write (77,9876) (it3(i),i=1,nmem)
    c9873            format(' INDXMAX:',/,(20i5))
    c9874            format(' DISTMAX:',/,(20f5.1))
    c9875            format(' IT2:',/,(20i5))
    c9876            format(' IT3:',/,(20i5))
    c                write (77,9877) (t1(i),i=1,nmem)
    c9877            format(' RAVINDX:',/,(20f5.1))
    c                write (77,*) 'ICNEXT=',icnext
                call mrgsortlist(icent(ic0+1),it1,it2,it3,it4,it5,nofcls)
    c       do ic=1,nofcls
    c         write (iout,7721) ic,ifirst(ic),ilast(ic),
    c    -      (iclst(ia),ia=ifirst(ic),ilast(ic))
    c7721     format(' ic=',i3,' IFIRST,ILAST=',2i5,/,(20i4))
    c       end do
    c*****Find all clusters in a network and sort atoms in a cluster by groups
    c     Input parameters:
    c     n0,n: Use atomnumbers (vertices) from n0 to n (inclusive)
    c     maxnode,maxgr: Array sizes - see dimension statement above
    c     Workspace arrays: it1,it2,it3
    c     Output parameters:
    c     nofcls: Number of clusters requested
    c     iclst,ifirst,ilast: The elements of the ig-th group are atoms
    c     (iclst(ia),ia=ifrst(ig),ilast(ig))
    c
    c     Description of the algorithm:
    c     k-means++
    c     Cluster centers are the COM's of the cluster members
    c     print *,'CLSTRS_K_COM start n0,n,nofcls=',n0,n,nofcls
    c     Initialization
          call indexit(index2d,1,MAX2D,0)
          call init_kmedoids(nofcls,index2d,ifirst,ilast,icent,rmsd2d,n0,n,
            call trnsfr(cent(1,i),c(1,icent(i)),3)
          call zeroiti(indxclst,n0-1,n)
    c       Find the nearest center to each node
    c       write (6,7711) 'Bef sort indxclst:',(indxclst(i),i=n0,n)
    c7711   format(' KMEANS ',a,(/i3,19i4))
    c         Find the new centers
              call trnsfr(cent_prev,cent,3*nofcls)
              call zeroit(cent,3*nofcls)
              call zeroiti(icent,0,nofcls)
                  cent(k,ic)=cent(k,ic)+c(k,i)
                    cent(k,ic)=cent(k,ic)/float(icent(ic))
                  call trnsfr(cent,cent_prev,3*nofcls)
          call indexit(it1,1,nnode,0)
    c    do i=n0,n
    c       w(i-n0+1)=indxclst(i)
    c     end do
          call mrgsrti(6,it1,indxclst(n0),nnode,it2,it3,it4,it5,n)
    c     write (6,7711) 'aft sort indxclst:',(it1(i),i=1,n-n0+1)
    c     write (6,7711) 'aft sort it1:',(it1(i),i=1,n-n0+1)
    c     write (6,1000) iter,(icent(k),k=1,nofcls)
    c      do ic=1,nofcls
    c       if (ifirst(ic) .le. ilast(ic)) then
    c        write (40,7721) ic,(iclst(ia),ia=ifirst(ic),ilast(ic))
    c7721     format(' ic=',i3,/,(20i4))
    c        else
    c          write (40,*) 'Cluster ',ic,' is empty'
    c      end do
    c1000  format(' Iteration ',i4,' Number of cluster members=',(10i5))
          common /rangen/ ixo
          character*31 qcent
    c     print *,'INIT_KMEDOIDS n0,n,mx2d=',n0,n,mx2d
          call askyn('Do you want to specify initial cluster centers',46,
              call getint(qcent,31,999999,1,n,icent(k),00)
              call askyn(
            call randpx(1,ran)
    c           Select cent(k) to have the largest rmin
    c#    MMC routine 462 lstmod: 07/15/02
    c*****Find all clusters in a network and sort atoms in a cluster by groups
    c     Input parameters:
    c     n0,n: Use atomnumbers (vertices) from n0 to n (inclusive)
    c     maxnode,maxneig,maxgr: Array sizes - see dimension statement above
    c     nneig(i) : number of neighbours (bonded) of node i
    c     ineig(j,i) : j-th neighbour of atom i
    c     Workspace array: idrop
    c     Output parameters:
    c     nofcls0: Number of disconnected clusters (groups) found previously
    c     nofcls: Number of disconnected clusters (groups) found
    c     iclst,ifirst,ilast: The elements of the ig-th group are atoms
    c     (iclst(ia),ia=ifrst(ig),ilast(ig))
    c
    c     Description of the algorithm:
    c     Pick the node with the largest # of neighbors (and lowest energy)
    c     Make it a cluster; remove it ; repeat
    c     Initialization
    c     print *,'CLSTRS_nn start n0,n,maxgr=',n0,n,maxgr
    c       Find node(s) with largest nn
    c       write (6,1004) nofcls,ifirst(nofcls),ilast(nofcls),nleft,nnmax,
    c    -    imin,nneig(imin)
    c       Now eliminate imin and its neighbour from ineig
            call trnsfi(idrop,ineig(1,imin),nneig(imin))
    c1004  format(' nofcls,ifirst(nofcls),ilast(nofcls),nleft,nnmax=',5i5,/,
    c     -  ' imin,nn(imin)=',2i5)
    c*****Find all clusters in a network and sort atoms in a cluster by groups
    c     Input parameters:
    c     n0,n: Use atomnumbers (vertices) from n0 to n (inclusive)
    c     idenclst=1: bond if share nnmin neighbors
    c     idenclst=2: bond is within cutoff and share at least nnmin neighbors
    c     idenclst=3: bond if within cutoff and each has at least nnmin neighbors
    c     Workspace arrays: it1,it2,it3,it4,it5,it6
    c     Output parameters
    c     nofcls: Number of clusters requested
    c     iclst,ifirst,ilast: The elements of the ig-th group are atoms
    c     (iclst(ia),ia=ifirst(ig),ilast(ig))
    c     Initialization
    c     print *,'CLSTRS_DENSITY nnmin=',nnmin
    c     Find neighbor number range
    c               Remove i-j bond
    c        do i=nfrst,n
    c          write (77,7943) i,(ing(j,i),j=1,nng(i))
    c7943      format(i5,' NN:',50i5)
    c        end do
            call clean_ng(nfrst,n,nng,ing,mx2d)
    c        print *,'CLEAN done'
    c        do i=nfrst,n
    c          write (78,7943) i,(ing(j,i),j=1,nng(i))
    c        end do
            call checknnlist(nfrst,n,ing,nng,nerr,mx2d)
            call clstrs(ing,nng,it1,nfrst,n,iclst,ifirst,ilast,0,nofcls,
    c       print *,'CLSTRS done'
    c       Remove bonds if the # of common neighbors is < nnmin
    c             Mark bond for removal by setting it to -ing(ii,i)
    c       Remove the negative and zero ing entries
            call clean_ng(nfrst,n,nng,ing,mx2d)
            call clstrs(ing,nng,it1,nfrst,n,iclst,ifirst,ilast,0,nofcls,
    c       Keep bond when # of common neighbors is >= nnmin; use nng2r list too
            call zeroiti(it3,nfrst-1,n)
    c       Move the R-2R neighbors together with the R neighbors
    c           Check for # of common neighbors
    c             Create bond
    c       Move the it3 bonds to the start in ing
            call trnsfi(nng(nfrst),it3(nfrst),n-nfrst+1)
            call clstrs(ing,nng,it1,nfrst,n,iclst,ifirst,ilast,0,nofcls,
    c     Sort clusters by the # of members
          call indexit(it1,1,nofcls,0)
          call mrgsrti(iout,it1,it6,nofcls,it2,it3,it4,it5,nofcls)
    c     print *,'MRGSRT done'
    c      write (6,8943) (it1(i),i=1,nofcls)
    c8943  format(' IT1:',/,(20i4))
    c      write (6,8944) (-t2(i),i=1,nofcls)
    c8944  format(' T2:',/,(20f4.0))
          call trnsfi(ifirst,it2,nofcls)
          call trnsfi(ilast,it3,nofcls)
          call trnsfi(iclst,it4,n-nfrst+1)
    c      do i=1,nofcls
    c        write (79,7891) i,ifirst(i),ilast(i),
    c     -    (iclst(j),j=ifirst(i),ilast(i))
    c7891    format(i5,' if,il=',i4,i5,' iclst:',50i4)
    c      end do
    c     Remove the eliminated nodes from the nn list of the rest
    c     Create a list of indices in index where the value of mask is isave
    c     For the atoms in c1(index1(i)), c2(index2(i)) get the best fit
    c     using the algorithm of Kabasch, Acta Cryst, A32, p922 (1976).
    c     c1 is assumed to be the reference structure.
    c     cc1,cc2 are used for temporary storage
    c     print *,'BESTOVERLAY LEVTEST,maxat,nat=',LEVTEST,maxat,nat
    c     Move both sets into COM frame
              cc1(k,i)=c1(k,index1(i))
              cc2(k,i)=c2(k,index2(i))
    c     write (77,7831) 'INDEX1',(index1(i),i=1,nat)
    c     write (77,7831) 'INDEX2',(index2(i),i=1,nat)
    c7831 format(1x,a,':',/,(20i5))
    c     write (77,7832) 'CC1',(i,(cc1(k,i),k=1,3),i=1,nat)
    c     write (77,7832) 'CC2',(i,(cc2(k,i),k=1,3),i=1,nat)
    c7832 format(1x,a,':',/,(i5,3f10.5))
            com1(k)=dcom1(k)/atwsum
            com2(k)=dcom2(k)/atwsum
              cc1(k,i)=cc1(k,i)-com1(k)
              cc2(k,i)=cc2(k,i)-com2(k)
    c     Calculate coordinate sums - Numerical recipes starts arrays from 1!!
    c     Find the eigenvectors a and eigenvalues mu of rr
          call dtred2(rr,3,3,diag,offdiag)
          call dtqli(diag,offdiag,3,3,rr,ierr)
    c     The rows of the matrix a are the eigenvectors
    c     Calculate rotation matrix (U in Kabasch's notation)
    c       Planar molecule: a(iz)=a(inz)xa(jnz), same for b
    c       Linear molecules
    c       First set the iz-th components to nonparalel to the inz-th
    c       a(jz)=a(inz)xa(iz)
    c       a(iz)=a(inz)xa(jz)
          call check_rotmat(rot,'KABSCH',6,ifail,LEVTEST)
    c#    MMC routine 257 lstmod: 03/20/01
    c*****Extract block averages from cumulative sum
    c#    MMC routine 255 lstmod: 05/02/13
    c*****Computes error bound with the method of batch means
          character*(*) label
          common /sortsat/ ixdat(MAXSORT),dat(MAXSORT),datsort(MAXSORT),
          character*12 uncorr,corr,low,decide
    c     Minimum critical values for correlation test:
    c     Maximum critical values for correlation test:
    c     print *,'BATCHMEAN START iout=',iout,' NPTS=',npts
          ci=999.0
    c     Calculate avg, sd on the full data set
    c     print *,'AVG,SD=',avg,sd,' N=',n
    c     Aggregate, if needed to max MAXSORT blocks
          call blankout(decide,1,12)
            call indexit(ixdat,1,n,0)
            call trnsfr(datsort,dat,n)
            call mrgsrt(iout,ixdat,datsort,n,ifst,ilst,itemp,temp,n)
    c           Use the table
    c         Use formula assuming normal distribution
              ci=sd/sqrt(float(n-1))    
    c*****Save the bits from ibitx into mapbit
    c       write (77,*)' write iw=',iw,' mapbi=',mapbi,
    c    -    ' ibdone,itodo=',ibdone,itodo
    c     if (ib .gt. 0) write (77,*)' writ iw=',iw,' mapbi=',mapbi
    c*****Extract the bits from mapbit into ibitx
    c       Set loop limits so that inside loops can run in parallel
    c       write (77,*)' read il=',il,' mapbi=',mapbi,
    c    -    ' ibdone,itodo=',ibdone,itodo
    c           if (abs(e(m))+dd.eq.dd) go to 2
    c           MM 10/21/2004
    c           if (iter .eq. 300) pause 'too many iterations'
                c=1.d0
                    c=g/f
                    c=c*s
                    c=1.d0/r
    c     print *,'--- Simplex optimization started'
    c     write (6,1711) y,((p(i,k),k=1,3),i=1,4)
    c1711  format(' Y=',4f10.5,/,(' P=',3f10.6))
            call euler(rot,p(ilo,1),p(ilo,2),p(ilo,3))
              call euler(rot,p(ilo,1),p(ilo,2),p(ilo,3))
    c       print *,'Reset ',iiter,' ilo=',ilo,' ihi=',ihi
    c     ypr=funk(pr)
    c       yprr=funk(prr)
    c       yprr=funk(prr)
    c         write (6,1711) y,((p(i,k),k=1,3),i=1,4)
    c             y(i)=funk(pr)
    c         Continue contracting if the newest worst point is higher than
    c         the old worst point
    c         print *,'yw=',yw,' ywrst=',ywrst
    c#    MMC routine 252 lstmod: 04/24/05
    c*****Prints the date and the time
          character*8 version
          character*(*) mark
          character*100 hostname
          character*12 today
          common /today_date/ ltoday,today
    c      : Intel Fortran code
    c     C@ AB : Absoft Fortran code
    c     C@ G7 : Gnu G77 Fortran code
    c     C@ UG : SGI IRIX Fortran code
    c     C@ HP : Hewlett-Packard Fortran code
    c     C@ AX : IBM AIX code
    c     C@ UX : Generic Unix code
    c     INTEL compiler: on Linux  only, requires -Vaxlib compilation option
          character*8 date_dat
          character*10 time_dat
          character*5 zone_dat
          character*5 months(12)
    c     Find out host name
          call blankout(hostname,1,100)
          call lastchar(hostname,lhostname,100)
          call zeroiti(idtvalue,0,8)
    c     Print time and date
            call date_and_time(date_dat,time_dat,zone_dat,idtvalue)
    c#    MMC routine 105 lstmod: 01/18/07
          character*4 lab
    c*****Check the integrity of stored constants
          call indexit(icntrl,1,20,0)
          close (10)
          close (10,status='delete')
          character*60 answers,tips
          common /helplist/ init,maxans,answers(1000),
          common /tiplist/ initt,maxtips,tips(100),
    c     The answer to question iq is in lined ifl(iq)-ill(ig);
    c     each line is linelen characters long
    c     print *,'EXPLANATION ihelp,itip=',ihelp,itip,' init,initt=',
    c    -  init,initt
    c       Initialize linelen, iff, ill
              call lastchar(answers(il),ilc,60)
    c       print *,'NHELP=',nhelp
    c       Initialize linelent, ifft, illt
              call lastchar(tips(il),ilc,60)
    c       print *,'NTIPS=',ntips
end module simulaid