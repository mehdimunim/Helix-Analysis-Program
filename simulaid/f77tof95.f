c*****Preprocessor to generate the source code of the program MMC
      character*80 line,line1,blankline,filename
      character*2 sym
      character*4 SYSTEM
      character*5 tfl(2)
      dimension id1(60),id2(30)
      common /values/ isize(70),iopt(30),maxo,maxs
c------------------------------------------------------------------------------
      do ic=1,80
        blankline(ic:ic)=' '
      end do
      write (6,2000)
c-----Input filename
200   write (6,2006)
      filename=blankline
      read (5,1000) filename
      icl=80
      do while (icl .gt. 1 .and. filename(icl:icl) .eq. ' ')
        icl=icl-1
      end do
      if (icl.eq. 1) go to 200
      namlen=icl
c-----Open files
      open(unit=10,status='old',file=filename(1:namlen),
     -   form='formatted',iostat=ios)
      if (ios .ne. 0) then
        print *,'Problem opening ',filename(1:namlen)
        go to 200
      end if
c     Add _f95 to the file name root
      ic=1
      do while (ic .lt. namlen .and. filename(ic:ic) .ne. '.')
        ic=ic+1
      end do
      if (ic .lt. namlen) filename(ic+4:namlen+4)=filename(ic:namlen)
      filename(ic:ic+3)='_f95'
      namlen=namlen+4
      open(unit=20,status='new',file=filename(1:namlen),
     -  form='formatted',iostat=ios)
      if (ios .ne. 0)
     -  open(unit=20,status='old',file=filename(1:namlen),
     -   form='formatted',iostat=ios)
      if (ios .ne. 0) then
        print *,'Problem opening ',filename(1:namlen)
        stop
      end if
      rewind 20
      nlong=0
      nline=0
      do while (.true.)
        line=blankline
        read (10,1000,end=999) line
        nline=nline+1
c       Find last nonblank
        icl=80
        do while (icl .gt. 1 .and. line(icl:icl) .eq. ' ')
          icl=icl-1
        end do
c       Change declaration formats
        if (line(1:16) .eq. '      character*') then
          if (line(17:19) .eq. '(*)') then
            line(23:icl+3)=line(20:icl)
            line(16:22)='(len=*)'
          else
            ic0=18
            do while (idigit(line(ic0:ic0)) .eq. 1)
              ic0=ic0+1
            end do
            icdel=ic0-17
            line(23+icdel:icl+5)=line(18+icdel:icl)
            line(20+1:20+icdel)=line(17:16+icdel)
            line(16:20)='(len='
            line(21+icdel:21+icdel+1)=') '
          end if
          icl=icl+5
        else if (line(1:12) .eq. '      real*8') then
          line(15:icl+1)=line(14:icl)
          line(11:14)='(8) '
          icl=icl+1
        else if (line(1:15) .eq. '      integer*2') then
          line(18:icl+1)=line(17:icl)
          line(14:17)='(2) '
          icl=icl+1
        end if
        write (20,1000) line(1:icl)
        if (icl .gt. 72 .and.
     -      line(1:1) .ne. 'c' .and. line(1:1) .ne. 'C') then
          nlong=nlong+1
          write (6,2001) line(1:icl)
        end if
      end do
999   print *
      print *,'Fortran 95 source code: ',filename(1:namlen)
      if (nlong .gt. 0) write (6,2001) nlong
      stop
1000  format(a)
2000  format(' Conversion from F77 to F95',10x,'(version: 10/25/05)')
2001  format(' Line is longer than 72 characters:',/,a)
2006  format(' Name of the source file=',$)
      end
      function idigit(in)
      character*1 in,idig(10)
      data idig /'0','1','2','3','4','5','6','7','8','9'/
      idigit=0
      do i=1,10
        if (in .eq. idig(i)) then
          idigit=1
          return
        end if
      end do
      return
      end
      subroutine askyn(q,lenq,iyn,idefans,ians)
      character*(*) q
      character*132 pline
      character*1 ans
      character*5 defans
c     idefans=-1: default no; idefans=+1: default yes
c     iyn=1: yes -> ians=1, no -> ians=0
c     iyn=0: yes -> ians=0, no -> ians=1
      if (idefans .eq. -1) then
        defans='[n] '
        lendef=4
      else if (idefans .eq. +1) then
        defans='[y] '
        lendef=4
      else
        defans=' '
        lendef=1
      end if
      pline(1:1)=' '
      pline(2:lenq+1)=q(1:lenq)
      pline(lenq+2:lenq+8)=' (y/n)'
      pline(lenq+9:lenq+8+lendef)=defans(1:lendef)
100   write (6,1000) pline(1:lenq+8+lendef)
      read (5,1001,end=99,err=99) ans
99    ians=0
      if (ans .ne. 'n' .and. ans .ne. 'N' .and. ans .ne. 'y' .and.
     -  ans .ne. 'Y') then
        if (idefans .eq. -1) then
          ans='n'
        else if (idefans .eq. 1) then
          ans='y'
        else
          print *,'Pls answer y or n'
          go to 100
        end if
      end if
      if (ans .eq. 'y' .or. ans .eq. 'Y') ians=1
      if (iyn .eq. 0) ians=1-ians
      return
1000  format(a,1x,$)
1001  format(a1)
      end
      subroutine getint(q,len,idefault,in)
      character*(*) q
      character*132 ansline,pline
      lenq=len
      pline(1:1)=' '
      pline(2:lenq+1)=q(1:lenq)
      if (idefault .ne. 999999) then
c       Put default on query
        pline(lenq+2:lenq+3)=' ['
        if (idefault .le. 9999) then
          write (pline(lenq+4:lenq+7),1002) idefault
          pline(lenq+8:lenq+9)=']='
          lenq=lenq+9
        else
          write (pline(lenq+4:lenq+11),1003) idefault
          pline(lenq+12:lenq+13)=']='
          lenq=lenq+13
        end if
      else
        pline(lenq+2:lenq+2)='='
        lenq=lenq+2
      end if
100   write (6,1000) pline(1:lenq)
      read (5,1001,end=99,err=99) ansline
      ii=1
      call nextchar(ansline,ii,132)
      i1=ii
      call nextblank(ansline,ii)
      i2=ii-1
      if (i1 .gt. i2) then
        if (idefault .eq. 99999) go to 100
        in=idefault
      else
        read (ansline(i1:i2),*,err=999) in
      end if
      return
99    in=idefault
      return
999   print *,'Invalid input for an integer'
      go to 100
1000  format(a,$)
1001  format(a132)
1002  format(i4)
1003  format(i8)
      end
      subroutine nextchar(line,ifc,len)
      character*(*) line
c     Finds the next nonblank in line
      if (ifc .gt. len-1) return
      ifc1=ifc
      do i=ifc1,len-1
        ifc=i
        if (line(i:i) .ne. ' ' .and. line(i:i) .ne. '  ') then
          if (line(i:i) .eq. '!') ifc=len
          return
        end if
      end do
      ifc=len
      return
      end
      subroutine nextblank(line,ifc)
      character*(*) line
c     Finds the next blank in line
      if (ifc .gt. 131) return
      ifc1=ifc
      do i=ifc1,131
        ifc=i
        if (line(i:i) .eq. ' ' .or. line(i:i) .eq. '   ') then
          if (line(i:i) .eq. '!') ifc=132
          return
        end if
      end do
      ifc=132
      return
      end
      block data
      common /values/ isize(70),iopt(30),maxo,maxs
      character*2 sizesym,optname
      character*10 sizename
      character*25 optlname
      character*38 sizelname
      common /names/ sizename(70),sizelname(70),sizesym(70),
     -  optname(30),optlname(30)
      data maxo /30/,maxs/70/
      data sizesym
     - /'MO','MA','SX','MM','UW','TN','VN','TE','VE','LS',
     -  'VW','ST','GR','TA','SV','VT','NA','TL','GT','GV',
     -  'DT','DM','RG','PG','WG','OR','GX','GY','GZ','CV',
     -  'W2','WS','WI','MI','TR','AT','UU','UV','TG','VG',
     -  'ND','DG','LG','GE','GQ','PP','PS','WM','TD','FE',
     -  'MH','LT','MD','DC','RC','MW','MS','NH','MG','HA',
     -  10*'**'/
      data optname/'NN','TN','NA','NL','TS','FR','DB','UX','UD','UG',
     -  'UR','PS','EF','  ','DM', 'HP','I2','VC','IB','G7','PG','FG',
     -  'RF','AB','  ','  ','  ','  ','ND','NV'/
      data sizename
     -/'maxmol    ','maxatmol  ','mxpxslt   ','maxsltmol ','maxwnnu   ',
     - 'maxnst    ','maxnsv    ','maxest    ','maxesv    ','maxloopslt',
     - 'maxwnnv   ','maxslt    ','maxgslt   ','maxtslt   ','maxslv    ',
     - 'maxss     ','maxat     ','maxtrgrgr ','maxstg    ','maxsvg    ',
     - 'maxsst    ','maxmst    ','maxgrid   ','maxpfgr   ','maxcggr   ',
     - 'maxorgr   ','maxxgr    ','maxygr    ','maxzgr    ','maxcav    ',
     - 'maxlin    ','maxausp   ','maxauit   ','maxavit   ','maxtors   ',
     - 'maxatyp   ','maxatypu  ','maxstmol  ','maxtgrid  ','maxwrgrid ',
     - 'maxgvv    ','maxdrgrid ','maxdagrid ','maxpegrid ','mxpxgslt  ',
     - 'maxcavps  ','maxpfsum  ','maxmatch  ','maxtagrid ','mxfeslt   ',
     - 'maxhunsite','mxlooptor ','mxdiffmol ','mxdiffcs  ','mxrescr   ',
     - 'maxwidslt ','maxphsmol ','maxhmneig ','maxmolfg  ','maxath    ',
     -  10*'          '/
      data optlname /
     -  'Solvent near-neighbor map','Solute near-neighbor map ',
     -  'Arithmetic bit-map code  ','Logical bit-map handling ',
     -  'Solute torque calculation','Force/torque calculations',
     -  'Debugging code           ','Generic Unix             ',
     -  'Berkeley Unix            ','SGI Unix                 ',
     -  'AIX Unix                 ','SGI auto parallelization ',
     -  'Intel Fortan calls       ','                         ',
     -  'MPI-Distributed memory   ','Hewlett-Packard          ',
     -  'Integer*2                ','Vectorized search        ',
     -  'Isobaric ensemble        ','Gnu Fortran77            ',
     -  'Cavity grid analysis     ','Field gradient calcs.    ',
     -  'Reaction-field correction','Absoft Fortran 90/95     ',
     -  4*'                         ',
     -  'Not MPI                  ','Non-vectorized search    '/
      data sizelname /
     -  'solvent molecules+1                   ',
     -  'atoms per solute molecule             ',
     -  'solute atoms for proximity analysis   ',
     -  'solute molecules                      ',
     -  'words for solute neighbour bit list   ',
     -  'nuclei on solute                      ',
     -  'nuclei on solvent                     ',
     -  'EPEN electrons on solute              ',
     -  'EPEN electrons on solvent             ',
     -  'number of solute molecules w loop move',
     -  'words for solvent neighbour bit list  ',
     -  'solute centers (all copies)           ',
     -  'solute groups (all copies)            ',
     -  'solute centers for torsion option     ',
     -  'solvent centers/solvent               ',
     -  'solute or solvent centers             ',
     -  'centers (atoms and pseudoatoms)       ',
     -  'solute groups within torsion groups   ',
     -  'solute centers with a general solvent ',
     -  'solvent centers in a general solvent  ',
     -  'solute centers for sensitivity analyss',
     -  'molecules for sensitivity analysis    ',
     -  'full g(r) and primary g(r) grid-points',
     -  'preferential sampling grid points     ',
     -  'coupling parameter distribution grids ',
     -  'energy difference distribution grids  ',
     -  'grids in the x dir for grid search    ',
     -  'grids in the y dir for grid search    ',
     -  'grids in the z dir for grid search    ',
     -  'cavities                              ',
     -  'adaptive US matching workspace        ',
     -  'stored probabilities                  ',
     -  'iterations allowed for adaptive US+1  ',
     -  'block average entries                 ',
     -  'torsions                              ',
     -  'atom types the program can store      ',
     -  'atom types in a given solute          ',
     -  'molecules or solute atoms             ',
     -  'total g(r) grid points                ',
     -  'grid points for solvent-solvent g(r)s ',
     -  'number of solvent-solvent g(r)s       ',
     -  'dipole correlation QCDF radial grids  ',
     -  'dipole correlation QCDF angular grids ',
     -  'solute-solvent PE QCDF energy grids   ',
     -  'different QCDFs                       ',
     -  'cavities with pref. sampl. weights    ',
     -  'preferential sampling weight sub sums ',
     -  'AUS iterations to match               ',
     -  'torsion angle distribution grids      ',
     -  'free energy solute atoms              ',
     -  'sites for Hungarian method matching   ',
     -  'torsion loops                         ',
     -  'molecules for diffusion and residence ',
     -  'structures for diffusion              ',
     -  'structures for residence time         ',
     -  'number of Widom solutes               ',
     -  'number of primary hydr shell molecules',
     -  'number of neighbors for full match try',
     -  'number of molecules for fg calculation',
     -  'number of atoms for site represent.   ',
     -  10*'                                      '/
c     This can be replaced by the data statements printed by MMC
c     using the PRCO SAVE command - it will not be overridden by the
c     statements in the driver
      data isize /70*2/
      data iopt /30*0/
      end
