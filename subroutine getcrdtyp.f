      subroutine getcrdtyp(inpt,filex,lenfile,ifiltyp,ixcluster,
     -  inoutlab,linoutlab)
      character*(*) filex
      character*5 crdext
      character*(*) inoutlab
      common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
     -  iommod,iommc,iommc4,iogro,iomol2,iomae,iocif,ioxxx,ioins,
     -  ionxyz,iosxyz,iosxyzrq,iograsp,iofull,lext(19),crdext(19)
      character*11 formatname
      common /formats/ iqdconv(20),formatname(19)
      character*1 ans
      character*5 ext
      character*132 line
c     print *,'GETCRDTYP filex=',filex(1:lenfile)
      ifiltyp=0
      lenext=lenfile
      if (lenext .ge. 5) then
c       Extract extension
        ic=lenext
        do while (filex(ic:ic) .ne. '.' .and. ic .gt. max0(1,lenext-6))
          ic=ic-1
        end do
        lenext=lenext-ic
        if (lenext .le. 5) ext(1:lenext)=filex(ic+1:ic+lenext)
      else
        ext(1:lenext)=filex(1:lenext)
      end if
      if (ext(1:3) .ne. 'dat' .and. lenext .le. 5) then
        do it=1,ioins
          if (crdext(it)(1:lenext) .eq. ext(1:lenext) .and.
     -        lext(it) .eq. lenext) ifiltyp=it
        end do
      end if
      if (ext(1:5) .eq. 'pdbqs') then
        ifiltyp=ioa3pdb
      else if (ext(1:5) .eq. 'pdbqt') then
        ifiltyp=ioa4pdb
      else if (lenext .eq. 3 .and. (ext(1:3) .eq. 'pdb' .or.
     -    ext(1:3) .eq. 'PDB' .or. ext(1:3) .eq. 'ent')) then
c       Specify PDB type
        write (6,2000) inoutlab(1:linoutlab),'PDB'
c       Check file for segid/chaind
        if (inpt .gt. 0) then
c         File exists, find out PDB version
          line(1:4)='    '
          do while (line(1:4) .ne. 'ATOM' .and. line(1:6) .ne. 'HETATM')
            read (inpt,1000,end=999) line
          end do
          rewind inpt
          ifiltyp=0
          if (line(22:22) .eq. ' ' .and. line(73:76) .ne. '    ') then
            write (6,2001) 'PDB','Charmm'
            ifiltyp=iocpdb
          else if (line(22:22) .ne. ' ' .and.
     -             line(73:76) .eq. '    ') then
            write (6,2001) 'PDB','Brookhaven'
            ifiltyp=iobpdb
          end if
          ians=0
          if (ifiltyp .gt. 0) call askyn('Is that OK',10,1,1,ians,00,0)
          if (ians .eq. 0) then
            call quiz(ans,ians,'b',' ',0,'PDB file type',13,0,5,6,0)
            if (ans .eq. 'b') ifiltyp=iobpdb
            if (ans .eq. 'c') ifiltyp=iocpdb
          end if
        else
c         New file - Brookhaven is assumed
          ifiltyp=iobpdb
          write (6,2003) 'Brookhaven'
        end if
      else if (ext(1:3) .eq. 'CRD') then
        write (6,2000) inoutlab(1:linoutlab),'Charmm CRD'
        if (inoutlab(1:linoutlab) .ne. 'output') then
c         Determine Charmm CRD type
          line (1:1)='*'
          do while (line(1:1) .eq. '*')
            read (inpt,1000,end=999) line
          end do
          read (inpt,1000) line
          call blankout(line,1,132)
          read (inpt,1000) line
          call lastchar(line,lc,132)
          rewind inpt
          if (lc .gt. 80) then
            write (6,2001) 'Charmm','extended'
            ifiltyp=iocha
          else
            write (6,2001) 'Charmm','standdard'
            ifiltyp=iochaex
          end if
          call askyn('Is that OK',10,1,1,ians,00,0)
        else
          ians=0
        end if
        if (ians .eq. 0) then
          call quiz(ans,ians,'o',' ',0,'Charmm CRD file type',20,0,5,6,
     -      0)
          if (ans .eq. 'o') ifiltyp=iocha
          if (ans .eq. 'e') ifiltyp=iochaex
        end if
      else if (ext(1:3) .eq. 'slt') then
c       Specify MMC version
        ifiltyp=iommc
        call askyn('Is the .slt file in the old format',34,1,-1,is4,134,
     -    0)
        if (is4 .eq. 1) ifiltyp=iommc4
      else if (ext(1:3) .eq. 'rep') then
        write (6,2000) inoutlab(1:linoutlab),'Macromodel/Xcluster'
        ifiltyp=iommod
        ixcluster=1
      end if
      if (ifiltyp .ne. 0) then
        write (6,2000) inoutlab(1:linoutlab),formatname(ifiltyp)
      else
100     call quiz(ans,ians,' ',inoutlab,linoutlab,'file format',11,0,5,
     -    6,0)
        ifiltyp=iqdconv(ians)
        if (ans .eq. 'x') then
c         Macromodel/Xcluster
          ifiltyp=iommod
          ixcluster=1
        end if
        if (ifiltyp .eq. 0) then
          print *,'PROGRAM ERROR: invalid answer character'
          go to 100
        end if
      end if
      return
999   write (6,2002) filex(1:max0(1,lenfile))
      stop
1000  format(a)
2000  format(' The ',a,' format is established as ',a)
2001  format(' The ',a,' format is found to be ',a)
2002  format(' ERROR: input file ',a,' does not have coordinate ',
     -  'records')
2003  format(' New file - ',a,' PDB selected')
      end
