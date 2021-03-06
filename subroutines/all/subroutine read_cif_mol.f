      subroutine read_cif_mol(nats,ian,iresno,isegno,nsegm,c,frocc,
     -  charge,altcol,inscol,index,line,lineread,title,ltitle,
     -  inptyp_cif,iofull,idcol,iruntyp,ialtcol,iinscol,ikeepfullalt,
     -  altnam,naltrec,naltdel,ninsres,ipredict,outfile,altfile,namleno,
     -  namlena,pdbid,ncol,asterisk,inpt,ierr,maxat)
      dimension ian(maxat),iresno(maxat),isegno(maxat),c(3,maxat),
     -  frocc(maxat),charge(maxat),index(maxat)
      character*1 altcol(maxat),inscol(maxat),altnam(50),asterisk
      character*132 line(maxat),linewr
      character*80 title
      character*200 linein,outfile,altfile
      character*4 pdbid
      character*2 iatnm2,an2
      common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
     -  mmatno(64),iatnm2(99)
      dimension limcol(2,30),lencol(30),ic_xyz(3)
      character*4 chain,chain_prev
      character*8 atomname
      character*10 rectyp
c     print *,'READ_CIF_MOL inptyp_cif=',inptyp_cif
      call setcol(inptyp_cif,ncol,idcol,ialtcol,iinscol,
     -  inamcol1,inamcol2,irescol1,irescol2,iccol1,iccol2,
     -  iresncol1,iresncol2,iseqncol1,iseqncol2,isegcol1,isegcol2,
     -  iresidcol1,iresidcol2,iqcol1,iqcol2,ipotcol1,ipotcol2,
     -  iocccol1,iocccol2,ichemcol1,ichemcol2,nrescol,nresncol,
     -  nsegcol,nnamcol,iofull)
      ierr=0
      call get_title_cif(inpt,linein,lineread,title,ltitle,pdbid)
      call askyn('Do you want to use the author (PDB) residue numbers',
     -  51,1,1,iusepdbres,000,0)
c     print *,'TITLE:',title(1:ltitle)
      do while (linein(1:11) .ne. '_atom_site.')
        read (inpt,1000,end=999) linein
      end do
      ic_seqno=0
      ic_atsym=0
      ic_atname=0
      ic_alt=0
      ic_ins=0
      ic_resname=0
      ic_segid=0
      call zeroiti(ic_xyz,0,3)
      ic_frocc=0
      ic_beta=0
      ic_qform=0
      ic_resno=0
      ic_resno_auth=0
      icol=1
      do while (linein(1:11) .eq. '_atom_site.')
        ic=1
        call nextblank(linein,ic,200)
        len=ic-1
        if (linein(12:len) .eq. 'id') then
          ic_seqno=icol
        else if (linein(12:len) .eq. 'type_symbol') then
          ic_atsym=icol
        else if (linein(12:len) .eq. 'label_atom_id') then
          ic_atname=icol
        else if (linein(12:len) .eq. 'label_alt_id') then
          ic_alt=icol
        else if (linein(12:len) .eq. 'pdbx_PDB_ins_code') then
          ic_ins=icol
        else if (linein(12:len) .eq. 'label_comp_id') then
          ic_resname=icol
        else if (linein(12:len) .eq. 'label_asym_id') then
          ic_segid=icol
        else if (linein(12:len) .eq. 'Cartn_x') then
          ic_xyz(1)=icol
        else if (linein(12:len) .eq. 'Cartn_y') then
          ic_xyz(2)=icol
        else if (linein(12:len) .eq. 'Cartn_z') then
          ic_xyz(3)=icol
        else if (linein(12:len) .eq. 'occupancy') then
          ic_frocc=icol
        else if (linein(12:len) .eq. 'B_iso_or_equiv') then
          ic_beta=icol
        else if (linein(12:len) .eq. 'pdbx_formal_charge') then
          ic_qform=icol
        else if (linein(12:len) .eq. 'label_seq_id') then
          ic_resno=icol
        else if (linein(12:len) .eq. 'auth_seq_id') then
          ic_resno_auth=icol
        end if
        icol=icol+1
        call blankout(linein,1,200)
        read (inpt,1000,end=999) linein
      end do
      if (ic_xyz(1)*ic_xyz(2)*ic_xyz(3) .eq. 0) then
        print *,'ERROR: missing coordinate record(s) is not found'
        ierr=1
      end if
      if (ic_resno .eq. 0) then
        print *,'ERROR: residue number record is not found'
        ierr=1
      end if
      if (ic_atname .eq. 0)
     -   print *,'WARNING: atom name record is not found'
      if (ic_resname .eq. 0)
     -   print *,'WARNING: residue name record is not found'
      if (ic_atsym .eq. 0)
     -   print *,'WARNING: atom symbol record is not found'
      lenerr=0
      chain_prev='    '
      iresno_prev=0
      iresno_curr=0
      iresno_auth_prev=0
      isg=0
      do while (linein(1:4) .eq. 'ATOM' .or. linein(1:6) .eq. 'HETATM')
        nats=nats+1
        index(nats)=nats
        call blankout(line(nats),1,132)
        do k=1,3
          c(k,nats)=999.9
        end do
        ian(nats)=-1
        iresno(nats)=-1
        occ=0.0
        tfac=0.0
        ic=1
        nc=0
        call lastchar(linein,lc,200)
        line(nats)(1:4)=linein(1:4)
        if (linein(1:6) .eq. 'HETATM') line(nats)(1:6)=linein(1:6)
        do while (ic .lt. lc)
          nc=nc+1
          call nextstring(linein,ic,limcol(1,nc),limcol(2,nc),200)
          lencol(nc)=limcol(2,nc)-limcol(1,nc)+1
        end do
        rectyp='XYZ coords'
        do k=1,3
          icr=ic_xyz(k)
          read(linein(limcol(1,icr):limcol(2,icr)),*,err=888,end=888)
     -      c(k,nats)
          lenccol=(iccol2-iccol1+1)/3
          iccol0=iccol1+(k-1)*lenccol-1
          call save_cif_rec(rectyp,10,'iiccol',6,limcol(1,icr),
     -      line(nats),linein,iccol0+1,iccol0+lenccol,lencol(icr),0,
     -      lenerr)
        end do
        if (ic_resno_auth .gt. 0) then
          rectyp='Res_aut_no'
          icr=ic_resno_auth
          read(linein(limcol(1,icr):limcol(2,icr)),*,err=888,end=888)
     -      iresno_auth
          if (iusepdbres .eq. 1) then
            iresno(nats)=iresno_auth
            call save_cif_rec(rectyp,10,'iiresncol',9,
     -        limcol(1,icr),line(nats),linein,iresncol1,iresncol2,
     -        lencol(icr),0,lenerr)
          end if
        end if
        if (ic_resno .gt. 0 .and. iusepdbres .eq. 0) then
          rectyp='Residue no'
          icr=ic_resno
          if (linein(limcol(1,icr):limcol(2,icr)) .eq. '.') then
c           Residue number has to be generated - temporarily set to zero
            iresno(nats)=-999999
          else
            read(linein(limcol(1,icr):limcol(2,icr)),*,err=888,end=888)
     -        iresno(nats)
            call save_cif_rec(rectyp,10,'iiresncol',9,
     -        limcol(1,icr),line(nats),linein,iresncol1,iresncol2,
     -        lencol(icr),0,lenerr)
            iresno_prev=iresno(nats)
          end if
        end if
        if (ic_atname .gt. 0) then
          icr=ic_atname
          call save_cif_rec('Atom name',9,'iinamcol',8,
     -      limcol(1,icr),line(nats),linein,inamcol1,inamcol2,
     -      lencol(icr),1,lenerr)
          atomname=line(nats)(inamcol1:inamcol2)
        end if
        if (ic_resname .gt. 0) then
          icr=ic_resname
          call save_cif_rec('Residue name',12,'iirescol',8,
     -      limcol(1,icr),line(nats),linein,irescol1,irescol2,
     -      lencol(icr),1,lenerr)
        end if
        if (ic_segid .gt. 0) then
          icr=ic_segid
          call save_cif_rec('Segment id',10,'iisegcol',8,
     -      limcol(1,icr),line(nats),linein,isegcol1,isegcol2,
     -      lencol(icr),1,lenerr)
          chain='    '
          chain(1:limcol(2,icr)-limcol(1,icr)+1)=
     -      linein(limcol(1,icr):limcol(2,icr))
          if (chain .ne. chain_prev) then
            isg=isg+1
            chain_prev=chain
          end if
          isegno(nats)=isg
        end if
        if (ic_atsym .gt. 0) then
          an2='  '
          icr=ic_atsym
          if (lencol(ic_atsym) .eq. 1) then
            an2(2:2)=linein(limcol(1,ic_atsym):limcol(2,ic_atsym))
          else
            an2(1:2)=linein(limcol(1,ic_atsym):limcol(2,ic_atsym))
            call uplow(an2(2:2),an2(2:2),1,noabc)
          end if
          call findname(an2,iatnm2,1,99,ian(nats),2)
          call save_cif_rec('Chemical sym',12,'iichemcol',9,
     -      limcol(1,icr),line(nats),linein,ichemcol1,ichemcol2,
     -      lencol(icr),0,lenerr)
        else if (ic_atname .gt. 0) then
          ian(nats)=ianum(atomname,1,lencol(ic_atname))
        end if
        if (ic_alt .gt. 0) then
          rectyp='Alt at id '
          icr=ic_alt
          if (linein(limcol(1,icr):limcol(2,icr)) .eq. '.' .or.
     -        linein(limcol(1,icr):limcol(2,icr)) .eq. '?')
     -      linein(limcol(1,icr):limcol(2,icr))=' '
          altcol(nats)=linein(limcol(1,icr):limcol(2,icr))
          call save_cif_rec(rectyp,9,'ialtcol',7,limcol(1,icr),
     -      line(nats),linein,ialtcol,ialtcol,lencol(icr),1,lenerr)
        end if
        if (ic_ins .gt. 0) then
          rectyp='Insrt code'
          icr=ic_ins
          if (linein(limcol(1,icr):limcol(2,icr)) .eq. '.' .or.
     -        linein(limcol(1,icr):limcol(2,icr)) .eq. '?')
     -      linein(limcol(1,icr):limcol(2,icr))=' '
          inscol(nats)=linein(limcol(1,icr):limcol(2,icr))
          if (inscol(nats) .ne. ' ') ninsres=ninsres+1
          call save_cif_rec(rectyp,9,'iinscol',7,limcol(1,icr),
     -      line(nats),linein,iinscol,iinscol,lencol(icr),1,lenerr)
c     subroutine save_cif_rec(label,llabel,clabel,lclabel,limcol,
c    -  line,linein,ic1,ic2,len,noadjust,ierr)
        end if
        if (ic_frocc .gt. 0) then
          rectyp='Occupancy '
          icr=ic_frocc
          read(linein(limcol(1,icr):limcol(2,icr)),*,err=888,end=888)
     -      frocc(nats)
          call save_cif_rec(rectyp,9,'iiocccol',8,limcol(1,icr),
     -      line(nats),linein,iocccol1,iocccol2,lencol(icr),0,lenerr)
        end if
        if (ic_beta .gt. 0) then
          rectyp='Beta      '
          icr=ic_beta
          read(linein(limcol(1,icr):limcol(2,icr)),*,err=888,end=888)
     -      charge(nats)
          call save_cif_rec(rectyp,4,'iiqcol',6,limcol(1,icr),
     -      line(nats),linein,iqcol1,iqcol2,lencol(icr),0,lenerr)
        end if
        if (ic_seqno .gt. 0) then
          rectyp='Sequence n'
          icr=ic_seqno
          read(linein(limcol(1,icr):limcol(2,icr)),*,err=888,end=888)
     -      ii
          call save_cif_rec(rectyp,10,'iiseqncol',9,limcol(1,icr),
     -      line(nats),linein,iseqncol1,iseqncol2,lencol(icr),0,lenerr)
        end if
        if (iresno(nats) .eq. -999999) then
c         Generate residue number
          if (iresno_auth .ne. iresno_auth_prev) then
c           New residue
            iresno_auth_prev=iresno_auth
            iresno_curr=iresno_prev+1
            iresno_prev=iresno_curr
          end if
          iresno(nats)=iresno_curr
        end if
        call altcolcheck(line(nats),idcol,iruntyp,ialtcol,
     -    ikeepfullalt,altnam,naltnam,naltrec,naltdel,ipredict,frocc,
     -    outfile,altfile,namleno,namlena,ncol,nats,asterisk,maxrec)
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
        read (inpt,1000,end=999) linein
      end do
999   if (lenerr .gt. 0) then
        print *,'Record-length errors were found'
        call askstop(0)
      end if
      call altdelcheck(nats,naltnam,naltrec,altnam,ipredict,iruntyp,
     -  altcol,ialtcol,idcol,naltdel,ncol,line,index,idrop,linewr,
     -  altfile,namlena,asterisk,maxrec)
      nsegm=isegno(nats)
c     do i=1,nats
c       write (78,7861) nats,naltrec,ialtcol,altcol(i),line(i)
c     end do
c7861 format(i6,' naltrec=',i4,' ialtcol=',i3,' altcol=',a,/,
c    -  ' line=',a,/)
      return
888   write (6,2002) rectyp,nats,linein(limcol(1,icr):limcol(2,icr)),
     -  (limcol(k,icr),k=1,2),icr
      stop
1000  format(a)
2002  format(' Invalid ',a,' entry for atom # ',i6,':',a,' lims:',2i4,
     -  ' icr=',i2)
      end
