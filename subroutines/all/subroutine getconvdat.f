      subroutine getconvdat(ifflist,ilflist,ncol,nconvdat,iconvtyp,
     -  n,line,index,inamcol1,inamcol2,nnamcol,irescol1,irescol2,
     -  igrpinfo,maxrec)
      character* 132 line(maxrec)
      dimension index(n)
      character*1 ans1
      character*80 inpline
      character*4 ires
      character*8 convdat
      character*200 sfilename
      common /savedat/ mxresdat,maxcondat,ifst(1000),ilst(1000),
     -  nres,iresgen,lsfilename,ires(1000),convdat(7,10000),sfilename
      character*8 atnam
      character*200 filenames(5)
      dimension iconvtyps(5),lfilenames(5)
      data lfilenames /18,15,11,13,13/,iconvtyps /2,1,1,1,2/
      data filenames /'cha_prot19_rtf.dat',
     -  'amb_prot_ua.dat','amb_all.dat','amb_all94.dat','cha_all22.dat'/
      character*57 RTFtyps
      common /explainRTF/ RTFtyps(4)
      if (ilflist .gt. 5) then
        print *,'PROGRAM ERROR: conversion list limit=',ilflist,' > 5'
        stop
      end if
      nleadsp=0
      nleaddigit=0
      igrpinfo=0
      do ia=1,n
        atnam(1:nnamcol)=line(index(ia))(inamcol1:inamcol2)
        if (atnam(1:1) .eq. ' ') nleadsp=nleadsp+1
        if (idigit(atnam(1:1),1) .eq. 1) nleaddigit=nleaddigit+1
      end do
      write (6,2006) RTFtyps
      write (6,2007) -1,'User-supplied',' ','RTF file'
      write (6,2007) 0,'User-supplied',
     -  ' conversion data file ','in Simulaid format'
      do i=ifflist,ilflist
        write (6,2007) i,'Simulaid-supplied ',
     -    'conversion data file ',filenames(i)(1:lfilenames(i))
      end do
120   call getint('The number corresponding to your choice',39,999999,
     -  0,ilflist,ifile,49)
      iconvtyp=0
      if (ifile .eq. 0) then
        call getname(sfilename,lsfilename,'Name of the conversion file',
     -    27,200,'',0,0,0,0)
        if (ncol .gt. 2) call quiz(ans1,iconvtyp,' ',' ',0,
     -    'PF label type',13,0,5,6,00)
      else if (ifile .ge. ifflist) then
        sfilename=filenames(ifile)
        lsfilename=lfilenames(ifile)
        iconvtyp=iconvtyps(ifile)
      else if (ifile .ne. -1) then
        print *,'Invalid number'
        go to 120
      end if
      nconvdat=0
      ntypedat=0
      if (ifile .ge. 0) then
c       Read in conversion information
        iuconv=80
        call openfile(iuconv,0,'conversion information',22,'old',
     -    sfilename,lsfilename,notfnd,1,1,0,0,0)
        write (6,2005) sfilename(1:lsfilename)
        do while (.true.)
          read (iuconv,1001,end=110) inpline
          nconvdat=nconvdat+1
c         write (77,*) nconvdat,inpline
          if (nconvdat .lt. maxcondat) then
            read (inpline,*,err=999) (convdat(k,nconvdat),k=1,ncol)
          else
            write (6,2001) maxcondat
            stop
          end if
        end do
999     write (6,2000) i,inpline
        stop
110     close (iuconv)
        print *,'Number of conversions read=',nconvdat
      end if
      if (ifile .ge. 0)
     -    call askyn('Do you want to read RTF file(s) (too)',
     -      37,1,-1,ichrtf,0,0)
      if (ichrtf .eq. 1 .or. ifile .lt. 0) then
        call getrtfdat(nconvdat,5,iconvtyp,ntypedat)
      end if
c      do i=1,nconvdat
c        write (77,7711) i,(convdat(k,i),k=1,ncol)
c7711    format(i5,5(2x,a8))
c      end do
c     Find the residue limits in convdat
      nres=1
      iresgen=0
      ires(nres)=convdat(1,1)(1:4)
      ifst(nres)=1
      do i=2,nconvdat
        if (convdat(1,i) .ne. convdat(1,i-1)) then
c         New residue found
          ilst(nres)=i-1
          nres=nres+1
          if (nres .gt. mxresdat) then
            write (6,2004) mxresdat,sfilename(1:lsfilename)
            stop
          end if
          ires(nres)=convdat(1,i)(1:4)
          ifst(nres)=i
          if (ires(nres) .eq. '    ') iresgen=nres
        end if
      end do
      ilst(nres)=nconvdat
c     write (77,*) 'iresgen,nconvdat,nres=',iresgen,nconvdat,nres
c     do ii=1,nres
c       write (77,*) ii,ires(ii),ifst(ii),ilst(ii)
c     end do
      nleadspcd=0
      nleaddigitcd=0
      do i=2,nconvdat
        if (convdat(2,1)(1:1) .eq. ' ') nleadspcd=nleadspcd+1
        if (idigit(convdat(2,1)(1:1),1) .eq. 1)
     -    nleaddigitcd=nleaddigitcd+1
        if (convdat(5,nconvdat)(1:5) .eq. 'GROUP') igrpinfo=1
      end do
      ipdbin=0
      if (float(nleaddigit)/float(n) .gt. 0.01) ipdbin=1
      ipdbcd=0
      if (float(nleaddigitcd)/float(nconvdat) .gt. 0.01) ipdbcd=1
      if (ipdbin .eq. 1 .and. ipdbcd .eq. 0) then
        write (6,2002)
        call askyn('Do you want to de-regularize them',33,1,1,idereg,0,
     -    0)
        if (idereg .eq. 1) then
          call fixrecform(line,index,n,3,inamcol1,inamcol2,
     -      irescol1,irescol2,maxrec)
        end if
      end if
      return
1001  format(a80)
2000  format(' ERROR: invalid conversion data record, record no:',
     -  i5,':',a)
2001  format(' ERROR: Data file has more than',i5,' lines')
2002  format(' Input names appear to be in regular PDB format ',
     -  '(e.g., 2HG1)',/,'but the conversion data do not (i.e., HG12)')
2004  format(' ERROR: Maximum number of different residues ',
     -  '(',i4,') is exceeded in file',/,5x,a)
2005  format(' Residue and atomname conversion table is in file ',a)
2006  format(/,' Conversion information can be obtained from ',/,
     -  3x,'a) conversion data files (see list below) and/or',/,
     -  3x,'b) RTF file(s) ',/,4(a,/),
     -  /,' Current options:')
2007  format(i3,2x,a,a,a)
      end
