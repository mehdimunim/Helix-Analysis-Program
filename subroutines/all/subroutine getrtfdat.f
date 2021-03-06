      subroutine getrtfdat(nconvdat,ncol,iconvtyp,ntypedat)
      character*80 linep,linech,lineprev
      character*4 ires,resnam,atnam,atnamo,potnam,attyp
      character*200 sfilename,rtffile,parmfile
      character*8 convdat
      common /savedat/ mxresdat,maxcondat,ifst(1000),ilst(1000),
     -  nres,iresgen,lsfilename,ires(1000),convdat(7,10000),sfilename
      character*1 lc
      dimension sig(1000),eps(1000),attyp(1000)
      data ictyp /0/,igrprev /0/,numat /0/
C     Read Charmm Amber or Gromacs topology (RTF) file
      nconvdat0=nconvdat
      iconvtyp=0
200   namlena=0
      call openfile(15,0,'RTF',3,'old',rtffile,namlena,notfnd,0,1,1,0,0)
c     First figure out if Amber, Charmm, or Gromacs
      ifound=0
      nlread=0
      do while (ifound .eq. 0)
        read (15,1000,end=999) linech
        ic=1
        if (linech(1:1) .eq. ' ') call nextchar(linech,ic,80)
        if (linech(ic:ic+3) .eq. 'DONE') then
          ifound=1
          ictyp=1
          write (6,2002) 'Amber'
        else if (linech(ic:ic+3) .eq. 'RESI' .or.
     -           linech(ic:ic+3) .eq. 'PRES') then
          ifound=1
          ictyp=2
          write (6,2002) 'Charmm'
        else if (linech(ic:ic) .eq. '[') then
          ifound=1
          ictyp=3
          write (6,2002) 'Gromacs'
          print *,'NOTE: only four-character names are used for now'
        end if
      end do
      if (iconvtyp .gt. 0) then
        if (iconvtyp .ne. ictyp) then
          write (6,2003)
          iconvtyp=0
        end if
      else
        iconvtyp=ictyp
      end if
      rewind 15
      if (iconvtyp .eq. 1) then
c       Read Amber RTF
        read (15,1000) linech
        read (15,1000) linech
        ic=1
        call nextblank(linech,ic,80)
        write (6,2004) linech(1:ic-1)
        read (15,1000) linech
        nlread=3
        do while (linech(1:4) .ne. 'STOP')
          read (15,1000,end=999) linech
          read (15,1000,end=999) linech
          nlread=nlread+2
          ic=1
          if (linech(1:1) .ne. ' ') icinc=-1
          if (linech(ic:ic) .eq. ' ') call nextchar(linech,ic,80)
          resnam=linech(ic:ic+3)
          read (15,1000) linech
          read (15,1000) linech
          read (15,1000) linech
          nlread=nlread+3
          atnam='xxxx'
          atnamo=atnam
c         write (77,*) 'resnam=',resnam
          do while (atnam .ne. '    ')
            call blankout(linech,1,80)
            read (15,1000) linech
            nlread=nlread+1
            ic=1
            call nextchar(linech,ic,80)
            if (ic .eq. 80) then
              atnam='    '
            else
              call nextblank(linech,ic,80)
              call nextchar(linech,ic,80)
              atnam=linech(ic:ic+3)
c             write (77,*) 'atnam=',atnam
              if (atnam .ne. atnamo .and. atnam .ne. 'DUMM') then
                nconvdat=nconvdat+1
                do icl=1,ncol
                  convdat(icl,nconvdat)='        '
                end do
c               Residue name
                convdat(1,nconvdat)(1:4)=resnam
c               Atom name
                convdat(2,nconvdat)(1:4)=atnam
c               Atom type
                call nextblank(linech,ic,80)
                call nextchar(linech,ic,80)
                convdat(3,nconvdat)(1:2)=linech(ic:ic+1)
                do is=1,8
                  call nextblank(linech,ic,80)
                  call nextchar(linech,ic,80)
                end do
                ic1=ic
                call nextblank(linech,ic,80)
c               Charge
                convdat(4,nconvdat)(1:ic-ic1)=linech(ic1:ic-1)
                write (convdat(5,nconvdat),2001) 1,1
                atnamo=atnam
c               write (77,7711)nconvdat,(convdat(ic,nconvdat),ic=1,ncol)
c7711            format(i5,7(1x,'|',a8,'|'))
              end if
            end if
          end do
c         Search for DONE
          do while (linech(1:4) .ne. 'DONE')
            read (15,1000) linech
            nlread=nlread+1
          end do
          read (15,1000,end=999) linech
          nlread=nlread+1
        end do
      else if (iconvtyp .eq. 2) then
c       Read Charmm RTF
        igr=0
        newgr=1
        do while (.true.)
          read (15,1000,end=999) linech
          nlread=nlread+1
          ic=1
          if (linech(1:1) .eq. ' ') call nextchar(linech,ic,80)
          if (linech(ic:ic+3) .eq. 'RESI' .or.
     -        linech(ic:ic+3) .eq. 'PRES') igr=0
          if (linech(ic:ic+3) .eq. 'GROU') then
            igr=igr+1
            newgr=1
          else if (linech(ic:ic+3) .eq. 'ATOM') then
            nconvdat=nconvdat+1
            do icol=1,ncol
              convdat(icol,nconvdat)='        '
            end do
c           Residue name
            convdat(1,nconvdat)(1:4)=resnam
c           Atom name
            ic=ic+4
            call nextchar(linech,ic,80)
            convdat(2,nconvdat)(1:4)=linech(ic:ic+3)
            ic=ic+4
c           Atom type
            call nextchar(linech,ic,80)
            convdat(3,nconvdat)(1:4)=linech(ic:ic+3)
            ic=ic+4
            ic0=ic
            lc=linech(ic:ic)
            do while (idigit(lc,2) .eq. 1 .and. ic .lt. ic0+14)
              ic=ic+1
              lc=linech(ic:ic)
            end do
            if (ic-ic0 .gt. 2) then
              call readreal(linech,ic0,ic-1,charge)
            else
              charge=0.0
              write (6,2006) (convdat(k,nconvdat)(1:4),k=1,3)
            end if
            write (convdat(4,nconvdat),2000) charge
            write (convdat(5,nconvdat),2001) newgr,igr
            if (newgr .eq. 1) newgr=0
          else if (linech(ic:ic+3) .eq. 'RESI' .or.
     -             linech(ic:ic+3) .eq. 'PRES') then
c           print *,'RES ',resnam,' done'
            resnam='    '
            ic=ic+4
            call nextchar(linech,ic,80)
            icl=ic
            call nextblank(linech,ic,80)
            resnam=linech(icl:ic-1)
          end if
        end do
      else
c       Read Gromacs RTF (.rtp) files
        lookres=1
        do while (.true.)
          call blankout(linech,1,80)
          read (15,1000,end=999) linech
          nlread=nlread+1
          ic=1
          if (linech(ic:ic) .eq. ' ') call nextchar(linech,ic,80)
100       if (lookres .eq. 1) then
            if (linech(ic:ic) .eq. '[') then
              call nextblank(linech,ic,80)
              call nextchar(linech,ic,80)
              if (linech(ic:ic+4) .eq. 'atoms') then
c               Residue found
                igrprev=0
                ic=1
c???            call getname4(ic,lineprev,resnam,80,1)
                call getname4(ic,lineprev,resnam,80,1)
                lookres=0
                numat=0
              end if
              lineprev=linech
            end if
          else
c           Read record
            if (linech(ic:ic) .eq. '[') then
c             Residue ended
              lookres=1
              lineprev=linech
              if (numat .eq. 1) go to 100
            else
              call getname4(ic,linech,atnam,80,1)
              call getname4(ic,linech,potnam,80,-1)
              if (ic .gt. 72) then
c               Residue ended
                lookres=1
              else
                nconvdat=nconvdat+1
                numat=numat+1
                do icol=1,ncol
                  convdat(icol,nconvdat)='        '
                end do
c               Residue name
                convdat(1,nconvdat)(1:4)=resnam
c               Atom name
                convdat(2,nconvdat)(1:4)=atnam
c               Atom type
                convdat(3,nconvdat)(1:4)=potnam
                call nextchar(linech,ic,80)
                icf=ic
                call nextblank(linech,ic,80)
                read (linech(icf:ic-1),*,err=900) charge
                write (convdat(4,nconvdat),2000) charge
                call nextchar(linech,ic,80)
c               print *,'names=',resnam,atnam,potnam
c               print *,'charge=',charge,' ic=',ic,' lic=',linech(ic:ic)
                read (linech(ic:ic),*,err=900) igr
                newgr=0
                if (igr .ne. igrprev) newgr=1
                write (convdat(5,nconvdat),2001) newgr,igr
                igrprev=igr
              end if
            end if
          end if
        end do
      end if
999   close (15)
      print *,'Finished reading',nlread,' lines'
      if (ictyp .eq. 0) then
        print *,'ERROR: file ',rtffile(1:namlena),' is not an RTF file'
      else
        print *,'Number of conversions read from this file=',
     -    nconvdat-nconvdat0
      end if
      call askyn('Do you have an other RTF file',29,1,0,morefile,0,0)
      if (morefile .eq. 1) go to 200
      if (ncol .gt. 5) then
        ntypedat=0
300     namlenp=0
        call openfile(16,0,'PARAM',5,'old',parmfile,namlenp,notfnd,
     -    0,1,1,0,0)
c       Extract L-J eps and sigma values for atomtypes
        if (ictyp .eq. 1) then
c         Amber
          print *,'Amber parameter file read is not implemented yet'
        else if (ictyp .eq. 2) then
c         Charmm
          linep(1:4)='    '
          do while (linep(1:4) .ne. 'NONB')
            read (16,1000) linep
          end do
          call lastchar(linep,ic,80)
c         Skip continuation line(s)
          lc=linep(ic:ic)
          do while (lc .eq. '-')
            read (16,1000) linep
            call lastchar(linep,ic,80)
            lc=linep(ic:ic)
          end do
c         Keep reading until blank line is found
          do while (ic .gt. 1)
            linep(1:1)=' '
            read (16,1000,end=998) linep
            icf=1
            call nextchar(linep,icf,80)
            if (icf .lt. 80) then
c             Read record
              ntypedat=ntypedat+1
              attyp(ntypedat)=linep(icf:icf+3)
c             Skip one number
              icf=icf+3
              call nextchar(linep,icf,80)
              call nextblank(linep,icf,80)
c             Read eps
              call nextchar(linep,icf,80)
              icl=icf
              call nextblank(linep,icl,80)
              call readreal(linep,icf,icl-1,epsilon)
c             Read sigma
              icf=icl
              call nextchar(linep,icf,80)
              icl=icf
              call nextblank(linep,icl,80)
              call readreal(linep,icf,icl-1,sigma)
              eps(ntypedat)=-epsilon
              sig(ntypedat)=2.0*sigma/2.0**(1.0/6.0)
c             write (77,*) 'ntypedat,eps,sig=',
c    -            ntypedat,eps(ntypedat),sig(ntypedat)
            end if
            call lastchar(linep,ic,80)
            if (ic .eq. 1 .and. linep(1:1) .eq. '!') ic=2
          end do
          print *,'Read ',ntypedat,' nonbonded parameters'
        else if (ictyp .eq. 3) then
c         Gromacs
          print *,'Gromacs parameter file read is not implemented yet'
        end if
998     close (16)
        call askyn('Do you have an other PARAM file',31,1,0,morefile,0,
     -    0)
        if (morefile .eq. 1) go to 300
c       write (77,2006) (attyp(i),eps(i),sig(i),i=1,ntypedat)
c2006  format(' Atom type   Eps     sigma',/,(5x,a4,2x,f8.4,f8.4))
        notype=0
        do ic=1,nconvdat
c         Find type, extract eps & sigma
          do it=1,ntypedat
            if (convdat(3,ic)(1:4) .eq. attyp(it)) then
              write (convdat(6,ic),2000) eps(it)
              write (convdat(7,ic),2000) sig(it)
              go to 980
            end if
            notype=notype+1
            write (convdat(6,ic),2000) 0.0
            write (convdat(7,ic),2000) 0.0
          end do
980       continue
        end do
        if (notype .gt. 0) print *,'WARNING: ',notype,
     -  ' atoms in the RTF file have no LJ parameters'
      end if
c     write (77,*) 'nconvdat=',nconvdat
c     do ncd=1,nconvdat
c       write (77,7711) ncd,(convdat(ic,ncd),ic=1,ncol)
c     end do
c7711  format(i5,7(1x,'|',a8,'|'))
      return
900   write (6,2005) linech,resnam,atnam,potnam,charge,igr
1000  format(a)
2000  format(f8.4)
2001  format('GROUP',i1,i2)
2002  format(' RTF file is in ',a,' format')
2003  format(' ERROR: RTF file  format is incompatible with the ',
     -  'conversion file read')
2004  format(' Amber database name:',a)
2005  format(' ERROR in Gromacs atom record :',/,a,/,' resnam=',
     -  a4,' atnam=',a4,' potnam=',a4,' charge=',f7.4,' igr=',i2)
2006  format(' WARNING: residue ',a,' atom ',a,' type ',a,
     -  ' has no charge entry (0 is assumed)')
      end
