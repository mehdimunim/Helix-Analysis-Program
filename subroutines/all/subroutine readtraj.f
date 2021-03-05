      subroutine readtraj(inpt,inptrajtyp,mmctrajtyp,n,naslv,nwatr,nslt,
     -  trajformatname,ntitltr,trtitle,trajnam,ltrajnam,natom,
     -  nfreeat,ifree,icntrl,c,ninconf,noutconf,increment,inpcrdtyp,
     -  ietotread,etot,ifail,ifirst,ilast,iconfsel,numsel,maxrepconf,
     -  nmc,lentest,tofac,maxconf,maxconfsel,maxrec)
      dimension c(3,maxrec),ifree(maxrec),icntrl(20),
     -  iconfsel(maxconfsel)
      character*11 trajformatname
      character*80 trtitle(32)
      character* 132 line
      character* 200 trajnam,trajnam_n(2)
      real*8 xtlabc,xtlabc0
      common /boxdat/ xtlabc(6),xtlabc0(6),box(3),box0(3),edge_gen(3,3),
     -  cell0(3,27),cell(3,27),cellalt(3,27),
     -  ncell,ioppbc,noboxinfoar,noboxinfow,noboxrep,
     -  istuner,iboxtypfound,ixcrd(3),ixang(3),ixyzhex(3),
     -  ixyzhextraj(3),isizewarn
      character*1 separatorchar
      common /filenuminfo/ iaskunderscore,separatorchar
      character*1 xyz
      common /axislab/ xyz(3)
      real*8 db(6)
      dimension istune(4),ltrajnam_n(2)
      data k /0/
c     print *,'READTRAJ inpt=',inpt,' nslt,n=',nslt,n,
c    -  ' inptrajtyp=',inptrajtyp
      ifail=0
      ncol=80
      etot=0.0
      ninconf=ninconf+1
      iverbconf=1
      if (ninconf .gt. maxrepconf) iverbconf=0
      ietotread=0
      if (inptrajtyp .eq. 1) then
        ncol=0
        if (ninconf .eq. 1) then
          limic11=0
          if (natom .gt. n) then
            n=natom
            write (6,1009) n
          end if
        else
          limic11=1
        end if
        if (icntrl(11) .gt. limic11) then
          read(inpt,err=1095,end=1099) (xtlabc(j),j=1,6)
          do k=1,3
            box(k)=xtlabc(ixcrd(k))*tofac
          end do
          noboxinfoar=0
        else
          noboxinfoar=1
        end if
        if (icntrl(9) .gt. 0 .and. ninconf .gt. 1) then
          do k=1,3
            read(inpt,err=1096,end=1099) (c(k,ifree(j)),j=1,nfreeat)
          end do
        else
          do k=1,3
            read(inpt,err=1098,end=1099) (c(k,j),j=1,natom)
          end do
        end if
        nmc=nmc+1
      else if (inptrajtyp .eq. 2) then
        ncol=0
        read (inpt,1012,err=1098,end=1099) ((c(k,i),k=1,3),i=1,n)
        if (noboxinfoar .eq. 0) then
          read (inpt,1012,err=1098,end=1099) box
          if (iboxtypfound .eq. 2) noboxinfoar=1
        end if
c       Just use input edge info
        do k=1,3
          xtlabc(ixcrd(k))=box(k)
          xtlabc(ixang(k))=90.0
        end do
        nmc=nmc+1
      else if (inptrajtyp .eq. 3) then
        nsvp=min0(naslv,3)
        ncol=0
        if (mmctrajtyp .eq. 1) then
c         Binary MMC
          ietotread=1
          read (inpt,err=2099,end=1099) nwatr,natr,db(1),nmc,nidmc,
     -      i3,i4,ia0,ia1,(db(i),i=2,6),r1
          etot=db(2)
          if (nwatr .lt. 0 .or. natr .lt. 0 .or. nmc .lt. 0
     -      .or. nwatr .gt. maxrec/naslv .or. natr .gt. maxrec) then
            write (6,1013)
            stop
          end if
          read (inpt,err=2099,end=1099)
     -      ((c(k,j),k=1,3),j=ia0,ia1),(((c(k,ia1+(i-1)*naslv+j),
     -      k=1,3),j=1,nsvp),i=1,nwatr)
          if (noboxinfoar .eq. 0) read (inpt,err=2099,end=1099) box
          if (istuner .gt. 0) then
c           Skip tuning info
            read (inpt,err=2099,end=1099) istune
            if (istune(1)+istune(2) .gt. 0)read (inpt,err=2099,end=1099)
            if (istune(3) .gt. 0) read (inpt,err=2099,end=1099)
            if (istune(4) .gt. 0) read (inpt,err=2099,end=1099)
          end if
          n=nslt+nwatr*naslv
        else if (mmctrajtyp .eq. 2 .or. mmctrajtyp .eq. 3) then
          read (inpt,1000,end=1099) line(1:80)
          read (line(1:80),1019,err=1097)
     -      natomp,nmolec,n,ia0,ia1,nmc,cplpar
          n=nslt+(nmolec-1)*naslv
          if (mmctrajtyp .eq. 2) then
c           ASCII MMC
            do i=ia0,ia1
              read (inpt,1000,end=1099) line(1:80)
              read (line(1:80),1017,err=1098) (c(k,i),k=1,3)
            end do
            do iw=2,nmolec
              do j=1,nsvp
                read (inpt,1000,end=1099) line(1:80)
                read (line(1:80),1017,err=1098)
     -            (c(k,nslt-2*naslv+iw*naslv+j),k=1,3)
              end do
            end do
          else if (mmctrajtyp .eq. 3) then
c           Annotated ASCII MMC
            do i=ia0,ia1
              read (inpt,1000,end=1099) line(1:80)
              read (line(1:80),1018,err=1098) (c(k,i),k=1,3)
            end do
            do iw=2,nmolec
              do j=1,nsvp
                read (inpt,1000,end=1099) line(1:80)
                read (line(1:80),1018,err=1098)
     -            (c(k,nslt-2*naslv+iw*naslv+j),k=1,3)
              end do
            end do
          end if
          read (inpt,1015,end=1099) db(1)
          read (inpt,1000,end=1099)
          etot=db(1)
          ietotread=1
          if (noboxinfoar .eq. 0) read (inpt,1010,end=1099) box
        else if (mmctrajtyp .eq. 4) then
c         PDB
          nmc=nmc+1
          nr=0
          line(1:3)='   '
          do while (line(1:3) .ne. 'END')
            call blankout(line,1,80)
            read (inpt,1000,end=1100) line(1:80)
            if (line(1:4) .eq. 'ATOM' .or. line(1:6) .eq. 'HETATM') then
              nr=nr+1
              read (line(31:54),1011,err=1098) (c(k,nr),k=1,3)
            else
              call checkforetot(6,line,ninconf,etot,ietotread,iverbconf)
            end if
          end do
1100      n=nr
          if (n .eq. 0) then
            if (nr .gt. 0) write (6,1020) ninconf
            ifail=1
            return
          end if
        else if (mmctrajtyp .eq. 5) then
c         Charmm CRD
          call blankout(line,1,80)
          read (inpt,1000,end=1099) line(1:80)
          do while (line(1:1) .eq. '*')
            call checkforetot(1,line,ninconf,etot,ietotread,iverbconf)
            read (inpt,1000,end=1099) line(1:80)
          end do
          read(line(1:5),1002,err=1098) n
          do i=1,n
            read (inpt,1000,end=1099) line(1:80)
            read (line(1:50),1003,err=1098) (c(k,i),k=1,3)
          end do
          nmc=nmc+1
        end if
      else if (inptrajtyp .eq. 6) then
c       Amber CDF

        nmc=nmc+1
      else if (inptrajtyp .eq. 4 .or.
     -  inptrajtyp .eq. 5 .and. ninconf .eq. 1) then
        ncol=132
        read (inpt,1000,err=1098,end=1099) trtitle(1)
        do i=1,n
          read (inpt,1001,err=1098,end=1099) line
          read (line(53:88),1016,err=1098) (c(k,i),k=1,3)
        end do
        nmc=nmc+1
      else if (inptrajtyp .eq. 5 .and. ninconf .gt. 1) then
        ncol=132
        read (inpt,1000,err=1098,end=1099) line
        do i=1,n
          read (inpt,1000,err=1098,end=1099) line
          read (line(6:41),1016,err=1098) (c(k,i),k=1,3)
        end do
        nmc=nmc+1
      end if
      if (inptrajtyp .gt. 1) then
c       Just in case, generate Charmm edges
        do k=1,3
          xtlabc(ixcrd(k))=box(k)*tofac
          xtlabc(ixang(k))=90.0
        end do
      end if
      return
1095  write (6,1004) ninconf,trajformatname,' xtal   ',' ',inpt
      if (ncol .gt. 0) write (6,1006) line(1:ncol)
      ifail=1
      return
1096  write (6,1004) ninconf,trajformatname,' coordinate ',xyz(k),inpt
      write (6,*) 'nfreeat,natom=',nfreeat,natom
      if (ncol .gt. 0) write (6,1006) line(1:ncol)
      ifail=1
      return
1097  write (6,1004) ninconf,trajformatname,' header ',' ',inpt
      if (ncol .gt. 0) write (6,1006) line(1:ncol)
      ifail=1
      return
1098  write (6,1004) ninconf,trajformatname,' coordinate ',xyz(k),inpt
      if (ncol .gt. 0) write (6,1006) line(1:ncol)
      ifail=1
      return
1099  if (ilast .eq. 999999) then
        write (6,1008) trajformatname,ninconf-1
        ilast=ninconf-1
c       print *,'n,nslt=',n,nslt
      else
c       Generate the next trajectory name
        call nextnames(trajnam,ltrajnam,trajnam_n,ltrajnam_n,ntraj)
        ntry=0
        do while (ntry .lt. ntraj)
          ntry=ntry+1
          call opentraj(c,0,inpt,inptrajtyp,n,ntitltr,trtitle,inpcrdtyp,
     -      ifirst,ilast,increment,maxconf,ninconf,noutconf,natom,
     -      nfreeat,ifree,icntrl,1,mmctrajtyp,trajnam_n(ntry),
     -      ltrajnam_n(ntry),'input trajectory',16,iconfsel,numsel,
     -      0,1,2,2,icellfound,notfnd,0,0,lentest,0,0,0,maxconfsel,
     -      maxrec)
          if (notfnd .eq. 0) then
            trajnam=trajnam_n(ntry)
            ltrajnam=ltrajnam_n(ntry)
            write (6,1014) trajnam(1:ltrajnam)
            ntry=ntraj
          end if
        end do
        if (notfnd .gt. 0) then
          write (6,1005) trajformatname,ninconf,
     -     (trajnam_n(ntry)(1:ltrajnam_n(ntry)),ntry=1,ntraj)
          if (ncol .gt. 0) write (6,1007) line(1:ncol)
          ifail=1
        end if
      end if
      return
2099  print *,'ERROR: invalid binary MMC record found - ending scan'
      ifail=1
      return
1000  format(a80)
1001  format(a132)
1002  format(i5)
1003  format(20x,3f10.0)
1004  format(' ERROR: configuration #',i5,' of the ',a6,
     -  ' trajectory has invalid',a,'input ',a,' inpt=',i2)
1005  format(' ERROR: ',a6,' trajectory unexpectedly ended',/,
     -  8x,'End was found while reading configuration #',i8,/,
     -  ' If the trajectory is in segments, possible names are ',/,
     -  (1x,a))
1006  format(' Line read:',/,a)
1007  format(' Last line read:',/,a)
1008  format(1x,a6,' trajectory was found to contain ',i5,
     -  ' configurations')
1009  format(' Number of atoms set to ',i6,' from the Charmm ',
     -  'trajectory')
1010  format(9x,3f10.0)
1011  format(3f8.0)
1012  format(10f8.3)
1013  format(' ERROR: invalid atom/molecule numbers indicate that',/,
     -  8x,'MMC binary trajectory does not have box info saved',/,
     -  ' - restart without asking for drawing the box')
1014  format(' Next trajectory segment ',a,' opened')
1015  format(19x,e16.6)
1016  format(3f12.5)
1017  format(3f15.5)
1018  format(5x,3f15.5)
1019  format(i6,9x,i6,8x,i6,7x,2i6,5x,i9,4x,f8.0)
1020  format(' ERROR: no ATOM or HETATM record was found when reading',
     -  ' configuration',i6)
      end
