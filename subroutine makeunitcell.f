      subroutine makeunitcell(inpt,inpfile,linpfile,inpcrdtyp,iotyp,
     -  iocif,c,ctemp,atw,n,iatnum,nresslt,nsegslt,segid4,atnames,
     -  resnames,iresno,froccin,charge,ixres,molsltlim,line,index,
     -  iasymbio,iuout,outfile,namleno,itemp, blankline,radtodeg,
     -  maxseg,maxrsd,maxat)
      character*(*) inpfile,outfile
      character*132 line(maxat),blankline
      character*4 segid4(maxseg)
      character*8 atnames(maxat),resnames(maxrsd)
      dimension c(3,maxat),ctemp(3,maxat),index(maxat),atw(maxat),
     -  iatnum(maxat),iresno(maxat),froccin(maxat),charge(maxat),
     -  ixres(maxat),molsltlim(3,maxseg),itemp(maxat)
      character*2 iatnm2
      common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
     -  mmatno(64),iatnm2(99)
      character*1 abc,digits,hexdigits
      common /charactersets/ ihex(25),abc(62),digits(14),hexdigits(25)
      character*1 xyz
      common /axislab/ xyz(3)
      dimension shiftasym(3,20),rotasym(3,3,20),abcabg(6),cmin(3),
     -  cmax(3),c0(3),edge(3),uxyz(3,3),cellxyz(3,3),vertex(3,8),
     -  crm_as(3,20),crmtest(3),ctest(3),shift(3),maskseg(20),
     -  iusedseg(20)
      character*1 segids(62),chainprev,chlist(20),chidused(60),chnew
      character*2 nm2
      character*132 linew
      character*200 cellfile,lineinp
      if (iasymbio .eq. 1) print *,'Generating a full unit cell'
      if (iasymbio .eq. 2) print *,'Generating biological oligomers'
      if (iasymbio .eq. 3)
     -  print *,'Generating all crsytal contacts of the asymetric unit'
      if ((iasymbio .eq. 1 .or. iasymbio .eq. 3) .and.
     -    inpcrdtyp .eq. iocif) then
        print *,
     -    'Sorry, PDBx/mmCIF has no full symmetry tracsformation info'
        return
      end if
      write (6,2000) n,outfile(1:namleno)
c      write (6,9782) 'A',(segid4(i),i=1,nsegslt)
c9782  format(1x,a,' SEGID4:',20(a4,'|'))
      call zeroit(abcabg,6)
      do i=1,62
        segids(i)=abc(i)
      end do
      do is=1,nsegslt
        do i=1,62
          if (segid4(is)(1:1) .eq. segids(i)) segids(i)=' '
        end do
      end do
      if (inpcrdtyp .eq. iocif .and. ispdb(iotyp) .eq. 1) then
c       Make atom names conform PDB convention
        nunknown=0
        do ia=1,n
          if (iatnum(ia) .gt. 0) then
            nm2=iatnm2(iatnum(ia))
            if (nm2(1:1) .eq. ' ') then
              if (atnames(ia)(1:1) .eq. nm2(1:1)) then
                do ic=8,2,-1
                  atnames(ia)(ic:ic)=atnames(ia)(ic-1:ic-1)
                end do
                atnames(ia)(1:1)=' '
              end if
            end if
          else
            nunknown=nunknown+1
          end if
        end do
        if (nunknown .gt. 0)
     -    print *,'NOTE: ',nunknown,' atoms had unknown atomic number'
      end if
      call askyn(
     -  'Is the unit cell information in the input STRUCTURE file',56,
     -  1,+1,isamef,000,0)
      if (isamef .eq. 1) then
        inptcell=inpt
        rewind inpt
        cellfile=inpfile
        len=linpfile
        if (ispdb(inpcrdtyp) .ne. 1 .and. inpcrdtyp .ne. iocif)
     -    write (6,2006) cellfile(1:len),'is not'
      else
        inptcell=98
        len=0
        call openfile(inptcell,0,'cell information',16,'old',cellfile,
     -    len,notfnd,0,1,1,0,0)
        if (cellfile(len-2:len) .ne. 'pdb' .and.
     -      cellfile(len-2:len) .ne. 'PDB' .and. ispdb(inpcrdtyp).eq. 1)
     -    write (6,2006) cellfile(1:len),'does not appear to be','PDB'
        if (cellfile(len-2:len) .ne. 'cif' .and. inpcrdtyp .ne. iocif)
     -    write (6,2006) cellfile(1:len),'does not appear to be','mmCIF'
      end if
      icellfound=-5
c     Find crytal info first
      do while (icellfound .le. 0)
        read (inptcell,1000,end=555) lineinp
        if (lineinp(1:6) .eq. 'CRYST1') then
          read (lineinp(7:54),*,err=888) abcabg
          icellfound=1
        else if (lineinp(1:6) .eq. '_cell.') then
          if (lineinp(7:14) .eq. 'length_a') then
            ic=15
            call nextstring(lineinp,ic,ic1,ic2,200)
            read (lineinp(ic1:ic2),*,err=888) abcabg(1)
            icellfound=icellfound+1
          else if (lineinp(7:14) .eq. 'length_b') then
            ic=15
            call nextstring(lineinp,ic,ic1,ic2,200)
            read (lineinp(ic1:ic2),*,err=888) abcabg(2)
            icellfound=icellfound+1
          else if (lineinp(7:14) .eq. 'length_c') then
            ic=15
            call nextstring(lineinp,ic,ic1,ic2,200)
            read (lineinp(ic1:ic2),*,err=888) abcabg(3)
            icellfound=icellfound+1
          else if (lineinp(7:17) .eq. 'angle_alpha') then
            ic=18
            call nextstring(lineinp,ic,ic1,ic2,200)
            read (lineinp(ic1:ic2),*,err=888) abcabg(4)
            icellfound=icellfound+1
          else if (lineinp(7:16) .eq. 'angle_beta') then
            ic=17
            call nextstring(lineinp,ic,ic1,ic2,200)
            read (lineinp(ic1:ic2),*,err=888) abcabg(5)
            icellfound=icellfound+1
          else if (lineinp(7:17) .eq. 'angle_gamma') then
            ic=18
            call nextstring(lineinp,ic,ic1,ic2,200)
            read (lineinp(ic1:ic2),*,err=888) abcabg(6)
            icellfound=icellfound+1
          end if
        end if
      end do
555   irect=0
      call unitmat(uxyz)
      if (icellfound .eq. 1) then
        write (6,2002) abcabg
        call trnsfr(edge,abcabg,3)
      else
        print *,'No unit cell information (CRST1/_cell. records) found'
        call getxyz('Edge length in the ',19,' direction (A)',14,
     -    999999.0,centinp,1,0)
        call getreal('ey-ez angle(in deg)',19,90.0,abcabg(4),1,0)
        call getreal('ex-ez angle(in deg)',19,90.0,abcabg(5),1,0)
        call getreal('ex-ey angle(in deg)',19,90.0,abcabg(6),1,0)
      end if
      if (abcabg(4) .eq. 90.0 .and. abcabg(5) .eq. 90.0 .and.
     -    abcabg(6) .eq. 90.0) irect=1
      if (irect .eq. 0) then
c       Generate cell's unit vectors
        cosxy=cos(abcabg(6)/radtodeg)
        uxyz(1,2)=cos(abcabg(6)/radtodeg)
        uxyz(2,2)=sin(abcabg(6)/radtodeg)
        uxyz(1,3)=cos(abcabg(5)/radtodeg)
        uxyz(2,3)=
     -    (uxyz(1,2)*uxyz(1,3)-cos(abcabg(4)/radtodeg))/uxyz(2,2)
        uxyz(3,3)=sqrt(1.0-uxyz(1,3)**2+uxyz(2,3)**2)
        write (6,2005) uxyz
      end if
      do i=1,3
        do k=1,3
          cellxyz(k,i)=edge(i)*uxyz(k,i)
        end do
      end do
      rewind inptcell
      nasym=0
      if (iasymbio .eq. 1 .or. iasymbio .eq. 3) then
        do while (.true.)
          read (inptcell,1000,end=999) lineinp
          if (lineinp(1:18) .eq. 'REMARK 290   SMTRY') then
            read (lineinp(19:19),*,err=888) k
            read (lineinp(20:23),*,err=888) nasym
            read (lineinp(24:53),*,err=888) (rotasym(k,i,nasym),i=1,3)
            read (lineinp(54:68),*,err=888) shiftasym(k,nasym)
          end if
        end do
999     if (nasym .eq. 0) then
          print *,'No transformation info was found in file ',
     -      cellfile(1:len)
          stop
        end if
        if (inptcell .ne. inpt) close (inptcell)
        call checkdim(n*nasym,maxat,'MAXREC',6,
     -    'Number of atoms in the extended system',38,0)
c       Now create the transformed coordinates
        print *,'Number of transformations found:',nasym
        do isym=2,nasym
          call rotate_c(c,n,rotasym(1,1,isym),c(1,(isym-1)*n+1),
     -      'SYMTRAN',7)
          call shiftmol(c(1,(isym-1)*n+1),n,shiftasym(1,isym),
     -      c(1,(isym-1)*n+1),1.0)
          do is=1,nsegslt
            molsltlim(1,(isym-1)*nsegslt+is)=(isym-1)*n+molsltlim(1,is)
            molsltlim(2,(isym-1)*nsegslt+is)=(isym-1)*n+molsltlim(2,is)
            if (molsltlim(3,is) .gt. 0) then
             molsltlim(3,(isym-1)*nsegslt+is)=(isym-1)*n+molsltlim(3,is)
            else
              molsltlim(3,(isym-1)*nsegslt+is)=molsltlim(3,is)
            end if
          end do
        end do
        ntot=n*nasym
c       Generate the unit cell vertices
        call zeroit(vertex,3*8)
        nv=0
        do ix=1,2
          do iy=1,2
            do iz=1,2
              nv=nv+1
              do k=1,3
                if (ix .eq. 2) vertex(k,nv)=vertex(k,nv)+cellxyz(k,1)
                if (iy .eq. 2) vertex(k,nv)=vertex(k,nv)+cellxyz(k,2)
                if (iz .eq. 2) vertex(k,nv)=vertex(k,nv)+cellxyz(k,3)
              end do
            end do
          end do
        end do
        nattot=nasym*n
        call compact_ucell(c,ctemp,itemp,molsltlim,ntot,nasym*nsegslt,
     -    vertex,cmin,cmax,c0,edge,uxyz,nshift)
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
        print *,'NSHIFT=',nshift,' NTOT=',ntot
        if (nshift .gt. 0) then
          do i=1,ntot
            write (line(index(i))(31:54),1001) (c(k,i),k=1,3)
          end do
          call extension(c,itemp,0,1,ntot,cmin,cmax,c0,1,0,vol)
        end if
        if (iasymbio .eq. 3) then
c         Gather all contacts
          print *,'Wait ... '
          do is=1,nasym
            call cofms(c(1,(is-1)*n+1),crm_as(1,is),n,atw)
          end do
          ncopy=1
          call trnsfr(ctemp,c,n*3)
          do is=2,nasym
            print *,'Testing asymetric unit #',is
            d_org=dist2(crm_as(1,1),crm_as(1,is))
            do ix=-2,2
              do iy=-2,2
                do iz=-2,2
                  call zeroit(shift,3)
                  do k=1,3
                    shift(k)=shift(k)+ix*cellxyz(k,1)
                    shift(k)=shift(k)+iy*cellxyz(k,2)
                    shift(k)=shift(k)+iz*cellxyz(k,3)
                  end do
c                 Test copy # is shifted by shift
                  call arrsum(crm_as(1,is),shift,crmtest,3)
c                 if (dist2(crm_as(1,1),crm_test) .lt. d_org*1.2) then
c                   See if this copy is in contact with asym #1
                    do ia=1,n
                      call arrsum(c(1,(is-1)*n+ia),shift,ctest,3)
                      do ja=1,n
                        if (dist2(ctest,c(1,ja)) .lt. 36.0) then
c                         Use this copy
                          ncopy=ncopy+1
                          call checkdim(n*ncopy,maxat,'MAXREC',6,
     -                      'Number of atoms in the extended system',38,
     -                      0)
                          call trnsfr(ctemp(1,(ncopy-1)*n+1),
     -                      c(1,(is-1)*n+1),n*3)
                          call shiftmol(ctemp(1,(ncopy-1)*n+1),n,shift,
     -                      ctemp(1,(ncopy-1)*n+1),1.0)
                          go to 100
                        end if
                      end do
                    end do
c                 end if
100               continue
                end do
              end do
            end do
          end do
          print *,'Number of contact copies=',ncopy-1
          call trnsfr(c,ctemp,3*n*ncopy)
        else
          ncopy=nasym
        end if
        ifs=1
        do is=nsegslt+1,ncopy*nsegslt
          do while (segids(ifs) .eq. ' ')
            ifs=ifs+1
          end do
          segid4(is)='    '
          segid4(is)(1:1)=segids(ifs)
          segids(ifs)=' '
        end do
        write (6,2001) ((segid4((isym-1)*nsegslt+is)(1:1),is=1,nsegslt),
     -    ' ',isym=1,ncopy)
        iseg=nsegslt
        ix_nslt=index(n)
        do i=1,ix_nslt
          call lastchar(line(i),lc,80)
          write (iuout,1000) line(i)(1:lc)
        end do
        do isym=2,ncopy
          chainprev=' '
          do i=1,n
            linew=line(index(i))(1:80)
            write (linew(31:54),1001) (c(k,(isym-1)*n+i),k=1,3)
            write (linew(07:11),1002) (isym-1)*n+i
            if (linew(22:22) .ne. chainprev) then
              write (iuout,1000) 'TER'
              iseg=iseg+1
              chainprev=linew(22:22)
            end if
            linew(22:22)=segid4(iseg)(1:1)
            call lastchar(linew,lc,80)
            write (iuout,1000) linew(1:lc)
          end do
        end do
        if (iasymbio .eq. 1) then
          icelldup=1
          do while (icelldup .gt. 0)
            call askyn(
     -        'Do you want to create an (other) duplicate cell',47,
     -        1,-1,icelldup,000,0)
            if (icelldup .eq. 1) then
              call zeroit(shift,3)
              do k=1,3
                call getint(
     -            'Shift factor (-1,0,1) in the '//xyz(k)//' direction',
     -            40,0,0,1,ixyzk,0)
                do l=1,3
                  shift(l)=shift(l)+ixyzk*cellxyz(l,k)
                end do
              end do
              nw=n*ncopy
              nsw=iseg
              call shiftmol(c,nw,shift,ctemp,1.0)
              do isym=1,ncopy
                write (iuout,1000) 'TER'
                chainprev=' '
                do i=1,n
                  linew=line(index(i))(1:80)
                  write (linew(31:54),1001)(ctemp(k,(isym-1)*n+i),k=1,3)
                  write (linew(07:11),1002) nw+(isym-1)*n+i
                  if (linew(22:22) .ne. chainprev) then
                    write (iuout,1000) 'TER'
                    iseg=iseg+1
                    chainprev=linew(22:22)
                  end if
                  linew(22:22)=segid4(iseg-nsw)(1:1)
                  call lastchar(linew,lc,80)
                  write (iuout,1000) linew(1:lc)
                end do
              end do
            end if
          end do
          call askyn('Do you want to write the cell vertices/edges too',
     -      48,1,+1,icellw,000,0)
          if (icellw .eq. 1) then
            print *,'Cell vertices will be chain X, all HE atoms'
            write (iuout,1000) 'TER'
            do nv=1,8
              write (iuout,2003) nasym*n+nv,nasym*nresslt+1,
     -          (vertex(k,nv),k=1,3)
            end do
            do inc=1,2
              incr=nattot+4*(inc-1)
              write (iuout,2004) incr+1,incr+2,incr+3
              write (iuout,2004) incr+2,incr+4
              write (iuout,2004) incr+3,incr+4
            end do
            do i=1,4
              write (iuout,2004) nattot+i,nattot+i+4
            end do
          end if
        end if
        write (iuout,1000) 'END'
      else
c       Generate biological oligomers
        nchread=0
        iread=0
        nsegtot=nsegslt
        nslttot=n
        do ic=1,nsegslt
          chidused(ic)=segid4(ic)(1:1)
        end do
        call zeroiti(iusedseg,0,nsegslt)
        if (ispdb(inpcrdtyp) .eq. 1) then
          write (iuout,2011) 'REMARK 350'
          nasym=0
          do while (nasym .eq. 0)
            if (iread .eq. 0) then
              call blankout(lineinp,1,200)
              read (inptcell,1000,end=777) lineinp
            else
              iread=0
            end if
            if (lineinp(1:18) .eq. 'REMARK 350 APPLY T') then
              call lastchar(lineinp,lc,200)
              ic=42
              call nextchar(lineinp,ic,200)
              nchread=0
              do while (ic .le. lc)
                nchread=nchread+1
                chlist(nchread)=lineinp(ic:ic)
                ic=ic+3
              end do
            else if (lineinp(1:18) .eq. 'REMARK 350   BIOMT') then
              k=0
              do while (lineinp(1:18) .eq. 'REMARK 350   BIOMT')
                read (lineinp(19:19),*,end=888,err=888) k
                read (lineinp(20:23),*,end=888,err=888) nasym
                read (lineinp(24:53),*,end=888,err=888)
     -           (rotasym(k,i,nasym),i=1,3)
                read (lineinp(54:68),*,end=888,err=888)
     -            shiftasym(k,nasym)
                call blankout(lineinp,1,200)
                read (inptcell,1000,end=777) lineinp
                iread=1
              end do
            end if
          end do
        else
c         Read CIF info
          write (iuout,2011) 'pdbx_struct_oper_list.'
          rewind inptcell
          do while (lineinp(1:38) .ne.
     -              '_pdbx_struct_assembly_gen.asym_id_list')
            call blankout(lineinp,1,200)
            read (inptcell,1000,end=777) lineinp
          end do
          call lastchar(lineinp,lc,200)
          ic=39
          call nextchar(lineinp,ic,200)
          nchread=0
          do while (ic .le. lc)
            nchread=nchread+1
            chlist(nchread)=lineinp(ic:ic)
            ic=ic+2
          end do
c         print *,'CHLIST=',(chlist(i),i=1,nchread)
          rewind inptcell
          do while (lineinp(1:23) .ne. '_pdbx_struct_oper_list.')
            call blankout(lineinp,1,200)
            read (inptcell,1000,end=777) lineinp
          end do
          nitems=0
          ic_isym=0
          ic_m11=0
          ic_m12=0
          ic_m13=0
          ic_m21=0
          ic_m22=0
          ic_m23=0
          ic_m31=0
          ic_m32=0
          ic_m33=0
          ic_v1=0
          ic_v2=0
          ic_v3=0
          call lastchar(lineinp,lc,200)
          do while (lineinp(1:23) .eq. '_pdbx_struct_oper_list.')
            nitems=nitems+1
            if (lineinp(24:25) .eq. 'id') then
              ic_isym=nitems
            else if (lineinp(24:35) .eq. 'matrix[1][1]') then
              ic_m11=nitems
            else if (lineinp(24:35) .eq. 'matrix[1][2]') then
              ic_m12=nitems
            else if (lineinp(24:35) .eq. 'matrix[1][3]') then
              ic_m13=nitems
            else if (lineinp(24:35) .eq. 'matrix[2][1]') then
              ic_m21=nitems
            else if (lineinp(24:35) .eq. 'matrix[2][2]') then
              ic_m22=nitems
            else if (lineinp(24:35) .eq. 'matrix[2][3]') then
              ic_m23=nitems
            else if (lineinp(24:35) .eq. 'matrix[3][1]') then
              ic_m31=nitems
            else if (lineinp(24:35) .eq. 'matrix[3][2]') then
              ic_m32=nitems
            else if (lineinp(24:35) .eq. 'matrix[3][3]') then
              ic_m33=nitems
            else if (lineinp(24:32) .eq. 'vector[1]') then
              ic_v1=nitems
            else if (lineinp(24:32) .eq. 'vector[2]') then
              ic_v2=nitems
            else if (lineinp(24:32) .eq. 'vector[3]') then
              ic_v3=nitems
            end if
            call blankout(lineinp,1,200)
            read (inptcell,1000,end=777) lineinp
          end do
          nasym=0
          do while (lineinp(1:1) .ne. '#')
            ic=1
            call lastchar(lineinp,lc,200)
            isym=1
            do i=1,nitems
              call nextstring(lineinp,ic,ic1,ic2,200)
c             print *,'I=',i,' IC1,2=',ic1,ic2
              if (i .eq. ic_isym) then
                read (lineinp(ic1:ic2),*,err=666) isym
              else if (i .eq. ic_m11) then
                read(lineinp(ic1:ic2),*,err=666) rotasym(1,1,isym)
              else if (i .eq. ic_m12) then
                read(lineinp(ic1:ic2),*,err=666) rotasym(1,2,isym)
              else if (i .eq. ic_m13) then
                read(lineinp(ic1:ic2),*,err=666) rotasym(1,3,isym)
              else if (i .eq. ic_m21) then
                read(lineinp(ic1:ic2),*,err=666) rotasym(2,1,isym)
              else if (i .eq. ic_m22) then
                read(lineinp(ic1:ic2),*,err=666) rotasym(2,2,isym)
              else if (i .eq. ic_m23) then
                read(lineinp(ic1:ic2),*,err=666) rotasym(2,3,isym)
              else if (i .eq. ic_m31) then
                read(lineinp(ic1:ic2),*,err=666) rotasym(3,1,isym)
              else if (i .eq. ic_m32) then
                read(lineinp(ic1:ic2),*,err=666) rotasym(3,2,isym)
              else if (i .eq. ic_m33) then
                read(lineinp(ic1:ic2),*,err=666) rotasym(3,3,isym)
              else if (i .eq. ic_v1) then
                read(lineinp(ic1:ic2),*,err=666) shiftasym(1,isym)
              else if (i .eq. ic_v2) then
                read(lineinp(ic1:ic2),*,err=666) shiftasym(2,isym)
              else if (i .eq. ic_v3) then
                read(lineinp(ic1:ic2),*,err=666) shiftasym(3,isym)
              end if
              if (ic2 .eq. lc .and. i .lt. nitems) then
                call blankout(lineinp,1,200)
                read (inptcell,1000,end=777) lineinp
                call lastchar(lineinp,lc,200)
                ic=1
              end if
            end do
            nasym=isym
            call blankout(lineinp,1,200)
            read (inptcell,1000,end=777) lineinp
          end do
          if (nasym .lt. 2) then
            write (6,2012) cellfile(1:len)
            return
          end if
        end if
        nats_transf=molsltlim(2,nsegslt)
        iad=0
        do isym=1,nasym
          write (6,2009) isym,((rotasym(i,j,isym),j=1,3),i=1,3),
     -      (shiftasym(i,isym),i=1,3)
c         Check if transformation causes actual change
          notransf=1
          do k=1,3
            if (shiftasym(k,isym) .ne. 0.0 .or.
     -          rotasym(k,k,isym) .ne. 1.0) notransf=0
            do kk=1,k-1
              if (rotasym(k,kk,isym) .ne. 0.0 .or.
     -            rotasym(kk,k,isym) .ne. 0.0) notransf=0
            end do
          end do
          if (nchread .eq. 0) then
c           No chain was specified - use all chains
            do ic=1,nsegslt
              maskseg(ic)=1
            end do
          else
c           Set mask for chains selected
            call zeroiti(maskseg,0,nsegslt)
            do ic=1,nchread
              do icc=1,nsegslt
                if (chlist(ic) .eq. segid4(icc)(1:1)) maskseg(icc)=1
              end do
            end do
          end if
          do ic=1,nsegslt
            if (maskseg(ic).eq. 1) then
              if (iusedseg(ic) .eq. 1 .and. notransf .eq. 1) then
                print *,'Chain ',segid4(ic)(1:1),' skipped'
                maskseg(ic)=0
              else
                iusedseg(ic)=1
              end if
            end if
          end do
c         print *,'MASKSEG=',(maskseg(i),i=1,nsegslt)
          do ic=1,nsegslt
            if (maskseg(ic) .eq. 1) then
              if (notransf .eq. 0) then
c               Generate new chain ID
                ich=0
                new=0
                do while (ich .lt. 62 .and. new .eq. 0)
                  new=1
                  ich=ich+1
                  do icc=1,nsegtot
                    if (abc(ich) .eq. chidused(icc)) new=0
                  end do
                end do
                chnew=abc(ich)
                chidused(nsegtot+1)=chnew
              else
                chnew=segid4(ic)(1:1)
              end if
              write (6,2008) segid4(ic)(1:1),chnew
c              write (6,9671) (chidused(i),i=1,nsegtot)
c9671          format(' CHIDUSED=',20a1)
c             Transform chain ic and write the corresponding PDB file
              is1=molsltlim(1,ic)
              is1o=is1
              icn=ic
              if (isym .gt. 1) then
                is1=nslttot+1
                icn=nsegtot+1
              end if
              natss=molsltlim(2,ic)-molsltlim(1,ic)+1
              isn=is1
              call rotate_c(c(1,is1o),natss,rotasym(1,1,isym),
     -          ctemp(1,isn),'BIOMOL',6)
              call shiftmol(ctemp(1,isn),natss,shiftasym(1,isym),
     -          ctemp(1,isn),1.0)
              if (isym .gt. 1) then
                molsltlim(1,icn)=nslttot+1
                molsltlim(2,icn)=nslttot+natss
                if (molsltlim(3,is) .gt. 0) then
                  molsltlim(3,icn)=nslttot+(molsltlim(3,ic)-is1+1)
                else
                  molsltlim(3,icn)=molsltlim(3,ic)
                end if
              end if
c             print *,'IC =',ic ,' if,l=',(molsltlim(k,ic),k=1,2)
c             print *,'ICN=',icn,' if,l=',(molsltlim(k,icn),k=1,2)
              if (inpcrdtyp .eq. iotyp) then
c               Just change the coordinates and chain ids
                do i=molsltlim(1,icn),molsltlim(2,icn)
                  io=i-molsltlim(1,icn)+molsltlim(1,ic)
                  linew=line(index(io))(1:80)
                  write (linew(31:54),1001)(ctemp(k,i),k=1,3)
                  linew(22:22)=chnew
                  call lastchar(linew,lc,80)
                  write (iuout,1000) linew(1:lc)
                end do
              else if (inpcrdtyp .eq. iocif) then
c               Create full new record
                do i=molsltlim(1,icn),molsltlim(2,icn)
                  call blankout(linew,1,132)
                  call createrec_fromcif(linew,iocif,iotyp,
     -              ctemp(1,i),iatnum(i-iad),atnames(i-iad),
     -              resnames(ixres(i-iad)),chnew,iresno(i-iad),
     -              froccin(i-iad),charge(i-iad),i,blankline)
                  call lastchar(linew,lc,132)
                  write (iuout,1000) linew(1:lc)
                end do
              end if
              if (ispdb(iotyp) .eq. 1) write (iuout,1000) 'TER'
              if (notransf .eq. 0) then
                nslttot=nslttot+natss
                nsegtot=nsegtot+1
              end if
            end if
          end do
          iad=iad+nats_transf
        end do
777     if (nasym .lt. 2) then
          write (6,2012) cellfile(1:len)
         return
        end if
        if (inptcell .ne. inpt) close (inptcell)
      end if
      return
888   print *,'Error in line:'
      print *,lineinp
      stop
666   write (6,2010) i,ic2,ic2,lineinp
      stop
1000  format(a)
1001  format(3f8.3)
1002  format(i5)
2000  format(' Transformations will be applied to ',i6,' atoms - ',
     -  'waters are ignored',/,
     -  ' PDB file generated will be written to ',a)
2001  format(' Chain IDs of the transformed segments: ',/,1x,(78a1))
2002  format(' Unit cell dimensions:',3f10.3,' A',/,
     -  ' Cell axis angles:    ',3f10.3,' deg')
2003  format('ATOM  ',i5,1x,' HE ',1x,'CEL',1x,'X',i4,1x,3x,3f8.3,
     -  '  1.0   0.0')
2004  format('CONECT',4i5)
2005  format(' Unit vectors of the non-rectangular unit cell:',/,
     -  ' ex=',3f9.5,/,' ey=',3f9.5,/,' ez=',3f9.5)
2006  format(' NOTE: cell information is in PDB syntax but the file',/,
     -  a,1x,a,' a ',a,' file')

2008  format(' Chain ',a,' transformed as chain ',a)
2009  format(' Transformation matrix #',i2,/,3(3f10.5,/),
     -  ' Shift vector:',3f8.3)
2010  format(' Invalid biological dimer data in line (i=',i2,' ic1,2=',
     -  2i4,'):',/,1x,a)
2011  format('REMARK Biological oligomers, specified by ',a,' records')
2012  format(' The ',a,' file did not have oligomer transformation ',
     -  'information',/,' - the input structure should be the right ',
     -  'oligomer')
      end
