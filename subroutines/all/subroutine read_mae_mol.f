      subroutine read_mae_mol(nats,nbonds,ian,iresn,isegno,nsegm,c,q,
     -  nn,in,index,line,inptyp_mae,iofull,inpt,iout,ierr,
     -  maxat,maxng)
      dimension ian(maxat),iresn(maxat),isegno(maxat),c(3,maxat),
     -  q(maxat),index(maxat)
      dimension nn(maxat),in(maxng,maxat)
      character*132 line(maxat)
      character*400 linein
      parameter (MAXCOL=50)
      character*20 items(MAXCOL)
      character*50 colid(MAXCOL)
      dimension litems(MAXCOL),lcolid(MAXCOL),iancount(99)
      character*2 iatnm2
      common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
     -  mmatno(64),iatnm2(99)
      character*4 segn4,segn4prev
c     print *,'READ_MAE inpt,iout,inptyp_mae=',inpt,iout,inptyp_mae
      call setcol(inptyp_mae,ncol,idcol,ialtcol,iinscol,
     -  inamcol1,inamcol2,irescol1,irescol2,iccol1,iccol2,
     -  iresncol1,iresncol2,iseqncol1,iseqncol2,isegcol1,isegcol2,
     -  iresidcol1,iresidcol2,iqcol1,iqcol2,ipotcol1,ipotcol2,
     -  iocccol1,iocccol2,ichemcol1,ichemcol2,nrescol,nresncol,
     -  nsegcol,nnamcol,iofull)
      iccol11=iccol1
      iccol12=iccol1+(iccol2-iccol1+1)/3
      iccol13=iccol12+(iccol2-iccol1+1)/3
      len=400
      ierr=1
      ix_x=0
      ix_y=0
      ix_z=0
      ix_rnum=0
      ix_rnam=0
      ix_anam=0
      ix_occ=0
      ix_tfac=0
      ix_atno=0
      ix_charge1=0
      ix_segn=0
      call read_head_mae('m_atom',6,nats,linein,len,inpt,iout)
c     print *,'NATS=',nats
c     Read item descriptors
      call read_colid_mae(ndcol,colid,lcolid,'atoms',5,linein,
     -  len,iend,inpt,iout,MAXCOL)
      if (iend .eq. 1) return
      do icol=1,ndcol
        lc=lcolid(icol)
        if (colid(icol)(1:lc) .eq. 'r_m_x_coord') then
          ix_x=icol
        else if (colid(icol)(1:lc) .eq. 'r_m_y_coord') then
          ix_y=icol
        else if (colid(icol)(1:lc) .eq. 'r_m_z_coord') then
          ix_z=icol
        else if (colid(icol)(1:lc) .eq. 'i_m_residue_number') then
          ix_rnum=icol
        else if (colid(icol)(1:lc) .eq. 's_m_pdb_residue_name' .or.
     -           colid(icol)(1:lc) .eq. 's_m_residue_name') then
          ix_rnam=icol
        else if (colid(icol)(1:lc) .eq. 's_m_chain_name') then
          ix_segn=icol
        else if (colid(icol)(1:lc) .eq. 's_m_pdb_atom_name' .or.
     -           colid(icol)(1:lc) .eq. 's_m_atom_name') then
          ix_anam=icol
        else if (colid(icol)(1:lc) .eq. 'r_m_pdb_occupancy') then
          ix_occ=icol
        else if (colid(icol)(1:lc) .eq. 'r_m_pdb_tfactor') then
          ix_tfac=icol
        else if (colid(icol)(1:lc) .eq. 'i_m_atomic_number') then
          ix_atno=icol
        else if (colid(icol)(1:lc) .eq. 'r_m_charge1') then
          ix_charge1=icol
        end if
      end do
c     print *,'ix_x,y,z=',ix_x,ix_y,ix_z
c     print *,'ix_rnam,rnum,anam,atno=',ix_rnam,ix_rnum,ix_anam,ix_atno
      if (ix_x*ix_y*ix_z .eq. 0)
     -   write (iout,*) 'WARNING: coordinate records are missing'
      if (ix_anam .eq. 0) write (iout,*)
     -   'WARNING: atom name record is missing - generic atom names',
     -  ' will be generated'
      do ia=1,nats
        do k=1,3
          c(k,ia)=999.9
        end do
        ian(ia)=-1
        iresn(ia)=-1
        occ=0.0
        tfac=0.0
        call blankout(line(ia),inamcol1,inamcol2)
      end do
      segn4prev='@#*&'
      nsegm=0
      nnoname=0
      nUNK=0
      do ia=1,nats
        call blankout(line(ia),1,ncol)
        call blankout(linein,1,len)
        read (inpt,1000,end=991) linein
        index(ia)=ia
        ic=1
        do icol=1,ndcol
          call nextstring(linein,ic,ic1,ic2,len)
          litems(icol)=ic2-ic1+1
          items(icol)(1:litems(icol))=linein(ic1:ic2)
c         write (6,*) 'IA,ICOL,IC1,IC2=',ia,icol,ic1,ic2
        end do
        icol=1
c        write (6,8791) ia,(ic,litems(ic),items(ic)(1:litems(ic)),ic=1,3)
c8791  format(' IA=',i2,3(' ic=',i2,' l=',i2,' item=',a))
        read (items(icol)(1:litems(icol)),*,end=992,err=992) ix
        if (ix .ne. ia) then
          write (iout,2002) 'atom',ia,ix
          return
        end if
        icol=ix_x
        if (items(max0(1,icol))(1:2) .ne. '<>' .and. icol .gt. 0) then
          read (items(icol)(1:litems(icol)),*,end=992,err=992) c(1,ia)
          line(ia)(iccol11:iccol11+litems(icol)-1)=
     -      items(icol)(1:litems(icol))
        end if
        icol=ix_y
        if (items(max0(1,icol))(1:2) .ne. '<>' .and. icol .gt. 0) then
          read (items(icol)(1:litems(icol)),*,end=992,err=992) c(2,ia)
          line(ia)(iccol12:iccol12+litems(icol)-1)=
     -      items(icol)(1:litems(icol))
        end if
        icol=ix_z
        if (items(max0(1,icol))(1:2) .ne. '<>' .and. icol .gt. 0) then
          read (items(icol)(1:litems(icol)),*,end=992,err=992) c(3,ia)
          line(ia)(iccol13:iccol13+litems(icol)-1)=
     -      items(icol)(1:litems(icol))
        end if
        icol=ix_rnum
        if (items(max0(1,icol))(1:2) .ne. '<>' .and. icol .gt. 0) then
          read (items(icol)(1:litems(icol)),*,end=992,err=992) iresn(ia)
          line(ia)(iresncol1:iresncol1+litems(icol)-1)=
     -      items(icol)(1:litems(icol))
        end if
        if (ix_rnam .eq. 0) then
          call blankout(line(ia),irescol1,irescol2)
          line(ia)(irescol1:irescol1+2)='LIG'
        else
          icol=ix_rnam
          if (items(max0(1,icol))(1:2) .ne. '<>' .and. icol .gt. 0)
c    -      resnam(ia)(1:litems(icol))=items(icol)(1:litems(icol))
     -      line(ia)(irescol1:irescol1+litems(icol)-1)=
     -        items(icol)(1:litems(icol))
          if (items(max0(1,icol))(1:3) .eq. 'UNK') nUNK=nUNK+1
        end if
        icol=ix_segn
        if (items(max0(1,icol))(1:2) .ne. '<>' .and. icol .gt. 0) then
          segn4='    '
          len4=min0(4,litems(icol))
          segn4(1:len4)=items(icol)(1:len4)
          if (segn4 .ne. segn4prev) then
            nsegm=nsegm+1
            segn4prev=segn4
          end if
          isegno(ia)=nsegm
        end if
        icol=ix_atno
        if (items(max0(1,icol))(1:2) .ne. '<>' .and. icol .gt. 0) then
          read (items(icol)(1:litems(icol)),*,end=992,err=992) ian(ia)
c         if (ian(ia) .le. 99)
c    -      line(ia)(ichemcol1:ichemcol1+1)=iatnm2(ian(ia))
c         print *,'IA=',ia,' INCHEMCOL1=',ichemcol1,' cnam=',
c    -      iatnm2(ian(ia))
        end if
        icol=ix_anam
        inamefound=0
        if (icol .gt. 0) then
          if (items(icol)(1:2) .ne. '<>')
     -      line(ia)(inamcol1:inamcol1+litems(icol)-1)=
     -        items(icol)(1:litems(icol))
          inamefound=1
        end if
        if (inamefound .eq. 0 .and. ian(ia) .gt. 0) then
          line(ia)(inamcol1:inamcol1+1)=iatnm2(ian(ia))
          inamefound=1
        end if
        if (inamefound .eq. 0) nnoname=nnoname+1
        icol=ix_occ
        if (items(max0(1,icol))(1:2) .ne. '<>' .and. icol .gt. 0) then
          read (items(icol)(1:litems(icol)),*,end=992,err=992) occ
          line(ia)(iocccol1:iocccol1+litems(icol)-1)=
     -      items(icol)(1:litems(icol))
        end if
        write (line(ia)(1:10),1001) ia
        icol=ix_charge1
        if (items(max0(1,icol))(1:2) .ne. '<>' .and. icol .gt. 0) then
          read (items(icol)(1:litems(icol)),*,end=992,err=992) q(ia)
          line(ia)(iqcol1:iqcol1+litems(icol)-1)=
     -      items(icol)(1:litems(icol))
        else
          icol=ix_tfac
          if (items(max0(1,icol))(1:2) .ne. '<>' .and. icol .gt. 0) then
            read (items(icol)(1:litems(icol)),*,end=992,err=992) q(ia)
            line(ia)(iqcol1:iqcol1+litems(icol)-1)=
     -        items(icol)(1:litems(icol))
          end if
        end if
      end do
      if (nnoname .gt. 0) then
        if (nnoname .lt. nats .and. nnoname .ne. nUNK) then
c         Just added hydrogen
          nha=0
          do ia=1,nats
            if (line(ia)(inamcol1:inamcol2) .eq. '    ') then
              nha=nha+1
              line(ia)(inamcol1:inamcol1)='H'
              ic=inamcol1+1
              call writeint(line(ia),ic,nha,nhlen)
            end if
          end do
        else
c         All (UNK) names are blank - generate them from atomic numbers
          call zeroiti(iancount,0,99)
          do ia=1,nats
            if (line(ia)(inamcol1:inamcol2) .eq. '    ') then
              iancount(ian(ia))=iancount(ian(ia))+1
              line(ia)(inamcol1:inamcol1+1)=iatnm2(ian(ia))
              if (iancount(ian(ia)) .le. 99) then
                ic=inamcol1+1
                call writeint(line(ia),ic,iancount(ian(ia)),nal)
              else
                write (iout,*) 'Atom count for ',iatnm2(ian(ia)),
     -            'exceeds 99'
              end if
            end if
          end do
        end if
      end if
c     print *,'NATS=',nats
      call read_head_mae('m_bond',6,nbonds,linein,len,inpt,iout)
c     print *,'NBONDS,maxng=',nbonds,maxng
      if (maxng .gt. 1) then
c       Read bond info
        call read_colid_mae(ndcol,colid,lcolid,'bonds',5,linein,
     -    len,iend,inpt,iout,MAXCOL)
        if (iend .eq. 1) return
        call zeroiti(nn,0,nats)
        call zeroiti(in,0,maxng*nats)
c       print *,'NDCOL=',ndcol
        do ib=1,nbonds
          call blankout(linein,1,len)
          read (inpt,1000,end=991) linein
          ic=1
          call nextstring(linein,ic,ic1,ic2,len)
          icol=1
          read (linein(ic1:ic2),*,end=992,err=992) ix
          if (ix .ne. ib) then
            write (iout,2002) 'atom',ib,ix
            if (iout .ne. 6) write (6,2002) 'atom',ib,ix
            ierr=1
            return
          end if
          ic=ic2+1
          call nextstring(linein,ic,ic1,ic2,len)
          icol=2
          read (linein(ic1:ic2),*,end=992,err=992) ib1
          ic=ic2+1
          call nextstring(linein,ic,ic1,ic2,len)
          icol=3
          read (linein(ic1:ic2),*,end=992,err=992) ib2
          nn(ib1)=nn(ib1)+1
          in(nn(ib1),ib1)=ib2
          nn(ib2)=nn(ib2)+1
          in(nn(ib2),ib2)=ib1
        end do
      end if
      ierr=0
      return
991   write (iout,2001) 'atom descriptors'
      if (iout .ne. 6) write (6,2001) 'atom descriptors'
      ierr=1
      return
992   write (iout,2000) colid(icol)(1:lcolid(icol)),
     -  items(icol)(1:litems(icol)),linein(1:50)
      if (iout .ne. 6) write (6,2000) colid(icol)(1:lcolid(icol)),
     -  items(icol)(1:litems(icol)),linein(1:50)
c     print *,'icol=',icol,' litems=',litems(icol)
      ierr=1
      return
1000  format(a)
1001  format(i10)
2000  format(' ERROR: invalid ',a,':',a,' in line',/,1x,a)
2001  format(' ERROR: run out of data while reading ',a)
2002  format(' ERROR: index misaligment when reading ',a,
     - ' i=',i6,' i(read)=',i6)
      end
