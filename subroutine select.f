      subroutine select(line,nrecdel,idcol,asterisk,n,nslt,index,ixres,
     -  isegcol1,isegcol2,iseqncol1,iseqncol2,inamcol1,inamcol2,
     -  irescol1,irescol2,iqcol1,iqcol2,charge,iatno,nneig,ineig,
     -  indexdel,iout,maxneig,maxrec)
c     Extract parts of the system
      dimension index(n),ixres(n),indexdel(n),charge(n),iatno(n),
     -  nneig(n),ineig(maxneig,n)
      character*1 asterisk,ans,lastans
      character*4 chainid,chainidn,atnam,atnamadj,resnam,resnamadj,
     -  namesel(200)
      character* 132 line(maxrec)
      data ans /' '/
c      write (6,7877) isegcol1,isegcol2,iseqncol1,iseqncol2,
c     -  inamcol1,inamcol2,irescol1,irescol2
c7877  format(' isegcol12=',2i4,' iseqncol12=',2i4,' inamcol12=',2i4,
c     -  ' irescol12=',2i4)
      nsegcol=isegcol2-isegcol1+1
      nrescol=irescol2-irescol1+1
      if (iseqncol2 .gt. iseqncol1) then
        call readint(line(index(nslt)),iseqncol1,iseqncol2,iseqno,1,1,
     -    irerr)
      else
        iseqno=nslt
      end if
      maxresno=ixres(n)
      write (6,1000) iseqno,maxresno
      if (iout .gt. 0 .and. iout .ne. 6)
     -  write (iout,1000) iseqno,maxresno
      nkeep=0
      call zeroiti(indexdel,0,n)
      do while (.true.)
        nrecdel0=nrecdel
        lastans=ans
        call quiz(ans,ians,' ',' ',0,'selecting option',16,0,5,6,0)
        if (ans .eq. 'c' .or. ans .eq. 's') then
          if (nsegcol .gt. 0) then
            call blankout(chainid,1,nsegcol)
            if (ans .eq. 'c') then
              call getname(chainid,len,'Chain ID to keep',16,4,'',0,0,0,
     -          0)
              if (iout .gt. 0)
     -          write (iout,2002) 'Keep only',chainid(1:nsegcol)
            else if (ans .eq. 's') then
              call getname(chainid,len,'Chain ID to drop',16,4,'',0,0,0,
     -          0)
              if (iout .gt. 0)
     -          write (iout,2002) 'Drop',chainid(1:nsegcol)
            end if
          else
            print *,'This input format does not have chain ID'
          end if
        else if (ans .eq. 'k') then
          nkeep=nkeep+1
          if (nkeep .gt. 1) then
            write (6,2003)
            call askstop(1)
          else
            call getrange(ifst,999999,ilst,999999,i,0,'atom to keep',
     -        12,n,0)
          end if
          if (iout .gt. 0) write (iout,2000) 'atoms',' ',ifst,ilst
        else if (ans .eq. 'd') then
          call getrange(ifst,999999,ilst,999999,i,0,'atom to drop',12,
     -      n,0)
          if (iout .gt. 0) write (iout,2001) 'atoms',' ',ifst,ilst
        else if (ans .eq. 'r') then
          nkeep=nkeep+1
          if (nkeep .gt. 1) then
            write (6,2003)
            call askstop(1)
          else
            call getrange(ifst,999999,ilst,999999,i,0,
     -        'residue to keep',15,maxresno,0)
          end if
          if (iout .gt. 0) write (iout,2000) 'residues',' ',ifst,ilst
        else if (ans .eq. 'e') then
          call getrange(ifst,999999,ilst,999999,i,0,
     -      'residue to drop',15,maxresno,0)
          if (iout .gt. 0) write (iout,2001) 'residues',' ',ifst,ilst
        else if (ans .eq. 't') then
c         Read atom name list to keep
          call getnamelist(namesel,4,nnamesel,'Atoms names to use',18,
     -      200)
          write (6,2004) 'Atom',(namesel(i),i=1,nnamesel)
        else if (ans .eq. 'u') then
c         Read residue name list to keep
          call getnamelist(namesel,nrescol,nnamesel,
     -      'Residue names to use',12,200)
          write (6,2004) 'Residue',(namesel(i),i=1,nnamesel)
        else if (ans .eq. 'v') then
c         Drop all solvents
          if (n .eq. nslt)
     -      print *,'There are no solvents in the system - check the ',
     -        'solvent residue name'
        else if (ans .eq. 'q') then
          if (nrecdel .eq. n) then
            print *,'NOTE: all atoms were deleted'
            call askstop(-1)
          else if (nrecdel .eq. 0 .and. lastans .ne. 'v') then
            print *,'NOTE: no atom was deleted'
            call askstop(-1)
          end if
          return
        else if (ans .eq. 'a') then
          if (iout .gt. 0) write (iout,2000) 'the alpha carbons'
        else if (ans .eq. 'b') then
          if (iout .gt. 0) write (iout,2000) 'the backbone'
        else if (ans .eq. 'h') then
          if (iout .gt. 0) then
            write (iout,2001) 'the aliphatic hydrogens'
            if (iqcol2 .gt. iqcol1) write (iout,*)
     -        'Hydrogen charges will be added to the carbon charge'
          end if
        else if (ans .eq. 'l') then
          if (iout .gt. 0) then
            write (iout,2001) 'all hydrogens'
            if (iqcol2 .gt. iqcol1) write (iout,*)
     -        'Hydrogen charges will be added to the heavy atom charge'
          end if
        end if
        ncfound=0
        ignoreseqno=0
        iseqask=1
        do ia=1,n
          idrop=0
          if (isegcol2 .ge. isegcol1) then
            chainidn=line(index(ia))(isegcol1:isegcol2)
          else
            chainidn='    '
          end if
          if (iseqncol2 .gt. iseqncol1 .and. ignoreseqno .eq. 0) then
            call readint(line(index(ia)),iseqncol1,iseqncol2,iseqno,1,1,
     -        irerr)
          else
            iseqno=ia
          end if
          if (iseqno .gt. 99990 .and. ignoreseqno .eq. 0) then
            ignoreseqno=1
            write (6,*) 'Sequence number read will be ignored'
          end if
          if (iseqno .ne. ia .and. ignoreseqno .eq. 0 .and.
     -        iseqask .eq. 1) then
            write (6,2005) iseqno,ia
            call askyn('Do you want to ignore the sequence number read',
     -        46,1,1,ignoreseqno,0,0)
            iseqask=0
          end if
          iresno=ixres(ia)
          atnam=line(index(ia))(inamcol1:inamcol2)
          call leftadjust4(atnam,atnamadj)
c          write (77,7778) ia,iseqno,iresno,atnam,atnamadj,
c     -      chainid(1:nsegcol),chainidn(1:nsegcol)
c7778      format(' ia,iseqno,iresno=',3i4,' atnam=',a,' atnamadj=',a,
c     -      ' chainid,n=',a,1x,a)
          if (ans .eq. 'c') then
            if (chainidn(1:nsegcol) .ne. chainid(1:nsegcol)) idrop=1
            if (idrop .eq. 0) ncfound=ncfound+1
          else if (ans .eq. 's') then
            if (chainidn(1:nsegcol) .eq. chainid(1:nsegcol)) idrop=1
            if (idrop .eq. 1) ncfound=ncfound+1
          else if (ans .eq. 'b') then
c           if (atnamadj(1:3) .ne. 'CA ' .and. atnamadj(1:2) .ne. 'C '
c    -          .and. atnamadj(1:2) .ne. 'O '
c    -          .and. atnamadj(1:2) .ne. 'N '
c    -          .and. atnamadj(1:2) .ne. 'H '
c    -          .and. atnamadj(1:3) .ne. 'HN ') idrop=1
             if (isbackbone(atnamadj(1:3),3) .eq. 0)  idrop=1
c            write (77,9781) ia,atnam,atnamadj(1:3),idrop
c9781        format(i6,' atnam=',a,' atnamadj=',a,' idrop=',i2)
          else if (ans .eq. 'a') then
            if (atnamadj(1:3) .ne. 'CA ') idrop=1
          else if (ans .eq. 't') then
            idrop=1
            do i=1,nnamesel
              if (atnamadj(1:4) .eq. namesel(i)) idrop=0
            end do
          else if (ans .eq. 'u') then
            resnam(1:nrescol)=line(index(ia))(irescol1:irescol2)
            call leftadjustn(resnam,resnamadj,nrescol)
            idrop=1
            do i=1,nnamesel
              if (resnamadj(1:nrescol) .eq. namesel(i)(1:nrescol))
     -          idrop=0
            end do
          else if (ans .eq. 'h' .or. ans .eq. 'l') then
            if (iatno(ia) .eq. 1) then
              if (nneig(ia) .gt. 0) then
                iac=ineig(1,ia)
                if (iatno(ineig(1,ia)) .eq. 6 .or. ans .eq. 'l') then
                  idrop=1
                  if (iqcol2 .gt. iqcol1)
     -              charge(iac)=charge(iac)+charge(ia)
                end if
              else
                print *,'Atom ',ia,' has no neighbour'
                call askyn('Do you want to keep it',22,1,-1,ikeep,0,0)
                if (ikeep .eq. 0) idrop=1
              end if
            end if
          else if (ans .eq. 'k') then
            if (iseqno .lt. ifst .or. iseqno .gt. ilst) idrop=1
          else if (ans .eq. 'd') then
            if (iseqno .ge. ifst .and. iseqno .le. ilst) idrop=1
          else if (ans .eq. 'r') then
            if (iresno .lt. ifst .or. iresno .gt. ilst) idrop=1
          else if (ans .eq. 'e') then
            if (iresno .ge. ifst .and. iresno .le. ilst) idrop=1
          else if (ans .eq. 'v') then
            if (ia .gt. nslt) idrop=1
          end if
c         if (idrop .gt. 0 .and.
c    -       line(index(ia))(idcol:idcol) .ne. asterisk) then
c           nrecdel=nrecdel+1
c           line(index(ia))(idcol:idcol)=asterisk
c         end if
          if (idrop .gt. 0 .and. indexdel(ia) .eq. 0) then
            nrecdel=nrecdel+1
            indexdel(ia)=1
            line(index(ia))(idcol:idcol)=asterisk
          end if
        end do
        if (ans .eq. 'c' .or. ans .eq. 's') then
          if (ncfound .gt. 0) print *,'Number of chain ',
     -      chainid(1:nsegcol),' atoms=',ncfound
          if (ncfound .eq. 0) print *,'NOTE: no chain ',
     -      chainid(1:nsegcol),' atoms were found'
        end if
        print *,(nrecdel-nrecdel0),' atoms deleted in this step'
        if (ans .eq. 'h' .and. iqcol2 .gt. iqcol1) then
c         Write back charges into line
          do ia=1,n
            call blankout(line(index(ia)),iqcol1,iqcol2)
            if (iqcol2-iqcol1 .gt. 6)  then
              write (line(index(ia))(iqcol2-8:iqcol2),1003) charge(ia)
            else
              write (line(index(ia))(iqcol2-5:iqcol2),1002) charge(ia)
            end if
          end do
        end if
      end do
1000  format(' Selecting atoms from the system',/,' Number of solute ',
     -  'atoms=',i6,' Largest solute residue number=',i6)
1002  format(f6.3)
1003  format(f9.5)
2000  format(' Keep only ',a,a,'in the range [',i6,',',i6,']')
2001  format(' Delete ',a,a,'in the range [',i6,',',i6,']')
2002  format(1x,a,' chain/segment ',a)
2003  format(' Only one range-keeping operation is meaningful.',/,
     -  ' For multiple ranges, use multiple deleting operationa')
2004  format(1x,a,' names to use:',/,
     -  (10(1x,a4)))
2005  format(' Atom index read (',i5,') differs from the sequence ',
     -  'number (',i8,')')
      end
