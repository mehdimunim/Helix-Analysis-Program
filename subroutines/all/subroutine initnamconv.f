      subroutine initnamconv(noconv)
c     Read in atom and residue name conversion tables, set up limits, etc
      character*200 filename,line
      character*8 namnam
      character*9 recnam
      common /convspec/ incon,ioutcon,ideoxy,namnam(10)
      character*3 ires0
      character*4 anamcontab,rnamcontab,atnam
      common /convdat/ nanamcon,nrnamcon,nres0,ifst(100),ilst(100),
     -  ncon,nconcol,anamcontab(1000,13),rnamcontab(13,100),ires0(100),
     -  icongro
      ideoxy=-1
      filename='pdb_nam.dat'
      namlen=11
      iuconv=70
      call openfile(iuconv,0,'namcon',6,'old',filename,namlen,
     -  notfnd,1,1,0,0,0)
      write (6,2001) filename(1:namlen)
c     Read in conversion information
      maxconcol=13
c     nconcol: number of colums containing conversion data
c     ncon: number of name conversions (PDB can have several columns)
      read (iuconv,*) ncon,nconcol
      if (nconcol .gt. maxconcol) then
        print *,'ERROR: Conversion data columns exceeds limit (',
     -    maxconcol,')'
        stop
      end if
c     namename: Names of the conventions
      read (iuconv,*) (namnam(i),i=1,ncon)
      do ic=1,ncon
        if (namnam(ic)(1:3) .eq. 'Gro') icongro=ic
      end do
      call chosename(namnam,ncon,'input',5,incon)
      call chosename(namnam,ncon,'output',6,ioutcon)
      if (incon .eq. ioutcon) then
        print *,'Names will be matched without conversion'
        noconv=1
        return
      else
        noconv=0
      end if
      if (namnam(ioutcon) .eq. 'Macromod') then
        print *,'Sorry, conversions TO Macromodel are not implemented'
        stop
      end if
      if (incon .eq. icongro .or. ioutcon .eq. icongro) write (6,2006)
      write (6,2004) namnam(incon),namnam(ioutcon)
      recnam='residname'
      do i=1,100
        read (iuconv,1000,end=992) line
        read (line,*,err=9921) (rnamcontab(j,i),j=1,nconcol)
c       read (iuconv,*,err=9921,end=992)
c       read (iuconv,*,end=992)
c    -    (rnamcontab(j,i),j=1,nconcol)
        if (rnamcontab(1,i) .eq. 'DONE') then
          nrnamcon=i
          go to 120
        end if
      end do
      go to 993
9921  write (6,2005) recnam,filename(1:namlen),i,line
      stop
120   irec=0
      recnam='atomnname'
      do i=1,100000
        read (iuconv,1000,end=991) line
        read (line,*,err=9921) (anamcontab(irec+1,j),j=1,nconcol)
c        write (77,2200) i,irec,(anamcontab(irec+1,j),j=1,nconcol)
c2200    format(' i,irec=',2i4,' anam=',12(1x,a4))
        if (anamcontab(irec+1,1) .eq. 'DONE') then
          nanamcon=irec
          go to 110
        end if
        atnam=anamcontab(irec+1,2)
        ndiff=0
        do j=3,nconcol
          if (atnam .ne. anamcontab(irec+1,j) .and.
     -        anamcontab(irec+1,j) .ne. '   ') ndiff=ndiff+1
        end do
c        write (77,7711) i,irec,nconcol,ndiff,anamcontab(irec+1,1),
c     -   anamcontab(irec+1,2)
c7711    format(' i,irec,nconcol,ndiff=',4i5,2a6)
        if (ndiff .gt. 0) irec=irec+1
        if (irec .eq. 1000) then
          print *,'Atom name table exceeds the maximum (1000)'
          stop
        end if
      end do
110   nres0=1
      ires0(nres0)=anamcontab(1,1)(1:3)
      ifst(nres0)=1
      do i=2,nanamcon
        if (anamcontab(i,1) .ne. anamcontab(i-1,1)) then
c         New residue found
          ilst(nres0)=i-1
          nres0=nres0+1
          ires0(nres0)=anamcontab(i,1)(1:3)
          ifst(nres0)=i
        end if
      end do
      ilst(nres0)=nanamcon
      close (iuconv)
c      do i=1,nres0
c        write (99,7123) i,ires0(i),ifst(i),ilst(i)
c7123    format(i3,1x,a3,3i5)
c      end do
c      do i=1,nanamcon
c        write (77,7755) i,(anamcontab(i,j),j=1,nconcol)
c7755    format(i5,10(2x,'*',a4,'*'))
c      end do
      return
991   print *,'ERROR: atom name table ended abruptly'
      stop
992   print *,'ERROR: residue name table ended abruptly'
      stop
993   print *,'ERROR: too many residues, increase size of ',
     -  'rnamcontab'
      stop
1000  format(a200)
2001  format(' Name-conversion file:',a)
2004  format(' Atom and residue name convention are changed from ',
     -  a,' to ',a)
2005  format(' ERROR: invalid ',a,' record in file ',a,
     -  ' record no=',i8,' :'/,a)
2006  format(' NOTE: only 4 characters will be used for conversion',/,
     -  7x,'Only protein conversions are complete for Gromacs')
      end
