      subroutine refineconv(n,cvffnam,iatnam,iatnum,iused,
     -  icvff,nneig,ineig,nhneig,nringnb,iringtyp,maxng,maxtyp,maxat)
      dimension nneig(maxat),ineig(maxng,n),iatnum(maxat),
     -  iatnam(maxat),icvff(maxat),nhneig(maxat),
     -  cvffnam(maxat),nringnb(maxat),iringtyp(maxat),iused(maxtyp)
      character*2 nnam
      character*4 cvffnam,iatnam,obl,ob,ct2,ct3
c     'Refine' conversion by differentiation of some cases
      data obl /'O   '/,ob /'OB  '/,ct2 /'CT2 '/,ct3 /'CT3 '/
c-----Determine ring membership type
c     print *,'REFINECONV n,maxng,maxtyp=',n,maxng,maxtyp
c     Count ring member neighbours
      do i=1,n
        iringtyp(i)=0
c       Set iringtyp(i) to -1 for ring atoms
        nnam=cvffnam(i)(1:2)
        if (nnam .eq. 'c5' .or. nnam .eq. 'cp' .or. nnam .eq. 'np'
     -      .or. nnam .eq. 'ni') iringtyp(i)=-1
        nringnb(i)=0
        nn=nneig(i)
        do j=1,nn
          nnam=cvffnam(ineig(j,i))(1:2)
          if (nnam .eq. 'c5' .or. nnam .eq. 'cp' .or. nnam .eq. 'np'
     -        .or. nnam .eq. 'ni') nringnb(i)=nringnb(i)+1
        end do
      end do
c     Set iringtyp(i) to 5 and 6, resp when ring size is known for sure
      do i=1,n
        if (cvffnam(i)(1:2) .eq. 'c5') iringtyp(i)=5
        if (cvffnam(i)(1:2) .eq. 'cp') iringtyp(i)=6
        if (cvffnam(i)(1:2) .eq. 'ni') iringtyp(i)=5
      end do
c     Set iringtyp(i) to 4 for ring junction atoms (to be refined later)
      do i=1,n
        if (nringnb(i) .gt. 2) iringtyp(i)=4
      end do
      nchange=0
      nrep=0
      do while (nchange .gt. 0 .or. nrep .eq. 0)
c     For ring nitrogens not set yet, deduce ring size from neighbours'
        nrep=nrep+1
        do i=1,n
          if (iatnum(i) .eq. 7 .and. iringtyp(i) .eq. -1) then
            nn=nneig(i)
            do in=1,nn
              if (iringtyp(ineig(in,i)) .eq. 5) then
                iringtyp(i)=5
                nchange=nchange+1
              end if
              if (iringtyp(ineig(in,i)) .eq. 6) then
                iringtyp(i)=6
                nchange=nchange+1
              end if
            end do
          end if
        end do
      end do
c     For ring junction, set type if possible
      do i=1,n
        if (iringtyp(i) .eq. 4) then
          nn=nneig(i)
          nn5=0
          nn6=0
          do in=1,nn
            if (iringtyp(ineig(in,i)) .eq. 5) nn5=nn5+1
            if (iringtyp(ineig(in,i)) .eq. 6) nn6=nn6+1
          end do
          if (nn5 .gt. 0 .and. nn6 .gt. 0) iringtyp(i)=1
          if (nn6 .eq. 0 .and. nn5 .eq. 3) iringtyp(i)=2
          if (nn5 .eq. 0 .and. nn6 .eq. 3) iringtyp(i)=3
        end if
      end do
c     Now -1: some ring atom; 0: No ring atom; 1:5/6 junction; 2: 5/5 junction
c     3: 6/6 junction; 4: some junction; 5:5-membered ring; 6:6-membered ring
      do i=1,n
        if (iatnam(i) .eq. obl) then
c         Separate carbonyl oxygen in acetic acid
          ic=ineig(1,i)
          if (nneig(i) .gt. 1 .or. iatnum(ic) .ne. 6) then
            write (6,100) nneig(i),iatnam(i),i
          else
            nnc=nneig(ic)
            if (nnc .eq. 2) then
              ij1=ineig(1,ic)
              ij2=ineig(2,ic)
              if (iatnum(ij1) .eq. 6 .and. iatnam(ij2) .eq. 'OH1 '
     -          .or. iatnum(ij2) .eq. 6 .and. iatnam(ij1) .eq. 'OH1 ')
     -          then
                iatnam(i)=ob
                iused(83)=iused(83)+1
                icvff(i)=83
              end if
            end if
          end if
        else if (cvffnam(i)(1:2) .eq. 'c ') then
          if (nhneig(i) .eq. 2) then
            iatnam(i)=ct2
            iused(84)=iused(84)+1
            icvff(i)=84
          else if (nhneig(i) .eq. 3) then
            iatnam(i)=ct3
            iused(85)=iused(85)+1
            icvff(i)=85
          end if
        else if (cvffnam(i)(1:2) .eq. 'cn' .and. nhneig(i) .eq. 2) then
          iatnam(i)=ct2
          iused(86)=iused(86)+1
          icvff(i)=86
        end if
      end do
      return
100   format(' ERROR: type O has',i2,' neighbours and the first is '
     -  ,a4,' i=',i4,/,' - check coordinate file')
      end
