      subroutine reseq(line,index,n,nslt,nsegm,iatfirst,
     -  iresfirst,iresidfirst,inptyp,iresno,resnames,numres,numslv,
     -  naslv,ninsres,resnamslv,iresnrestart,iresidrestart,nconfig,
     -  ireseq,ireseqdef,maxrsd,maxrec)
      dimension iresno(maxrec)
      character* 132 line(maxrec)
      character*8 resnam,resnames(maxrsd),resnamslv
      character*4 segnam
      dimension index(n)
      character*5 crdext
      common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
     -  iommod,iommc,iommc4,iogro,iomol2,iomae,iocif,ioxxx,ioins,
     -  ionxyz,iosxyz,iosxyzrq,iograsp,iofull,lext(19),crdext(19)
      common /solvres/ iaskdiffres,idiffres
      common /logging/ logfile,ipredict
c     Make sure residue numbers and atom numbers are consecutive
      call setcol(inptyp,ncol,idcol,ialtcol,iinscol,
     -  inamcol1,inamcol2,irescol1,irescol2,iccol1,iccol2,
     -  iresncol1,iresncol2,iseqncol1,iseqncol2,isegcol1,isegcol2,
     -  iresidcol1,iresidcol2,iqcol1,iqcol2,ipotcol1,ipotcol2,
     -  iocccol1,iocccol2,ichemcol1,ichemcol2,nrescol,nresncol,
     -  nsegcol,nnamcol,iofull)
      nrescol=irescol2-irescol1+1
      nresncol=iresncol2-iresncol1+1
      nsegcol=isegcol2-isegcol1+1
c     print *,'RESEQ nslt,naslv,nconfig=',nslt,naslv,nconfig
c     print *,'reseq  n=',n,' nrescol,nresncol,nsegcol=',
c    -  nrescol,nresncol,nsegcol
c     print *,irescol1,irescol2,iresncol1,iresncol2,
c    -  isegcol1,isegcol2
c     print *,' iseqncol1,2=',iseqncol1,iseqncol2
c     print *,'iresidcol1,iresidcol2=',iresidcol1,iresidcol2
      if (ireseq .eq. -1) return
      if (nconfig .eq. 1) then
        if (ninsres .gt. 0) write (6,1101) ninsres
        call askyn('Do you want to adjust atom and residue numbers',
     -    46,1,ireseqdef,ireseq,0,0)
        if (ireseq .eq. 1) then
          call getint('Atom number of the first atom',29,
     -      1,1,0,iatfirst,48)
          call getint('Residue number of the first residue',35,
     -      1,1,0,iresfirst,48)
c         print *,'First atom and residue numbers=',iatfirst,iresfirst
          if (ischarmm(ioutyp) .eq. 1) then
            call getint('Residue id of the first residue',31,
     -        iresfirst,1,0,iresidfirst,48)
            print *,'First residue id=',iresidfirst
          else
            iresidfirst=0
          end if
          if (nsegm .gt. 1 .or. ipredict .eq. 1) then
c           Allow for consecutive residue numbers over segments
            call askyn(
     -       'Do you want to restart residue numbering at each segment',
     -        56,1,-1,iresnrestart,0,0)
            if (ischarmm(inptyp) .eq. 1) then
            call askyn(
     -      'Do you want to restart residue id numbers at each segment',
     -        57,1,+1,iresidrestart,0,0)
            else
              iresidrestart=0
            end if
            if (iresnrestart+iresidrestart .gt. 0)
     -        print *,'Residue ranges for each segment are listed above'
          end if
        end if
      end if
      if (ireseq .ne. 1) return
      numslv=0
      resnam='        '
      resnam(1:nrescol)=line(index(1))(irescol1:irescol2)
      iresnoprev=iresno(1)
      resnames(1)=resnam
      if (nsegcol .gt. 0)
     -  segnam(1:nsegcol)=line(index(1))(isegcol1:isegcol2)
c     print *,'resnam,resnum,segnam,n=',resnam,resnum,segnam,n
      numres=iresfirst
      iiresnum=iresfirst
      iiresid=iresidfirst
      ncha=0
c     print *,'C nsegcol=',nsegcol
      do ia=1,n
        line(index(ia))(iinscol:iinscol)=' '
        if (line(index(ia))(irescol1:irescol2) .ne. resnam(1:nrescol)
     -     .or. iresno(ia) .ne. iresnoprev) then
          numres=numres+1
          iiresnum=iiresnum+1
          iiresid=iiresid+1
          resnam(1:nrescol)=line(index(ia))(irescol1:irescol2)
          iresnoprev=iresno(ia)
          resnames(numres)=resnam
          if (resnamslv(1:nrescol) .eq. resnam(1:nrescol))
     -      numslv=numslv+1
        else if (iaskdiffres .eq. 0 .and. nconfig .eq. 1) then
          if (ia .eq. nslt+naslv+1) then
c           Second solvent is the same residue as the first
            call askyn(
     -        'Do you want the solvents to be different residues',49,
     -        1,1,idiffres,0,0)
            iaskdiffres=1
          end if
        end if
        if (ia .gt. nslt+naslv .and. idiffres .eq. 1) then
          if (mod(ia-1-nslt,naslv) .eq. 0) then
            numres=numres+1
            iiresnum=iiresnum+1
            iiresid=iiresid+1
            numslv=numslv+1
          end if
        end if
        if (nsegcol .gt. 0) then
          if (segnam(1:nsegcol) .ne. line(index(ia))(isegcol1:isegcol2))
     -                           then
c           New segment was found
            segnam(1:nsegcol)=line(index(ia))(isegcol1:isegcol2)
            if (iresnrestart+iresidrestart .gt. 0)
     -        write (6,*) 'Segment ',segnam(1:nsegcol)
            if (iresnrestart .eq. 1) then
              call getint('Residue number of the first residue',35,
     -          iresfirst,1,0,iresfirst,48)
              iiresnum=iresfirst
            end if
            if (iresidrestart .eq. 1) then
              call getint('Residue id of the first residue',31,
     -          iresidfirst,1,0,iresidfirst,48)
              iiresid=iresidfirst
            end if
          end if
        end if
c       Sequence number (if used)
        if (iseqncol1 .gt. 0) then
          call readint(line(index(ia)),iseqncol1,iseqncol2,iaorg,1,1,
     -      irerr)
          ianew=ia+iatfirst-1
          if (iaorg .ne. ianew) ncha=ncha+1
          if (inptyp .eq. iochaex) then
            write (line(index(ia))(iseqncol1:iseqncol2),1001) ianew
          else
            write (line(index(ia))(iseqncol1:iseqncol2),1000) ianew
          end if
        end if
        if (numres .gt. maxrsd) then
          write (6,2000) maxrsd,numres
          stop
        end if
c       Check for change in residue number/id
        call readint(line(index(ia)),iresncol1,iresncol2,iiresorg,2,1,
     -    irerr)
        if (iiresorg .ne. iiresnum) ncha=ncha+1
        if (ischarmm(inptyp) .eq. 1) then
c         Residue id
          call readint(line(index(ia)),iresidcol1,iresidcol2,iiresorg,2,
     -      1,irerr)
          if (iiresorg .ne. iiresid) ncha=ncha+1
        end if
        if (inptyp .eq. iocha ) then
          write (line(index(ia))(iresncol1:iresncol2),1000)
     -      mod(iiresnum,100000)
          if (iiresnum .eq. 100000)
     -       write (6,1100) 'residue number',100000
          write (line(index(ia))(iresidcol1:iresidcol2),1006)
     -      mod(iiresid,100000)
          if (iiresid .eq. 100000)
     -       write (6,1100) 'residue id',100000
        else if (inptyp .eq. iochaex) then
          write (line(index(ia))(iresncol1:iresncol2),1001) iiresnum
          write (line(index(ia))(iresidcol1:iresidcol2),1002) iiresid
        else if (inptyp .eq. iommc ) then
          write (line(index(ia))(iresncol1:iresncol2),1000)
     -      mod(iiresnum,100000)
          if (iiresnum .eq. 100000)
     -     write (6,1100) 'residue number',100000
        else
          write (line(index(ia))(iresncol1:iresncol2),1006)
     -      mod(iiresnum,10000)
          if (iiresnum .eq. 10000) write (6,1100) 'residue number',10000
        end if
      end do
c     Add '3' to 'TIP'
c     print *,'numres=',numres
      do ir=1,numres
        if (resnames(ir)(1:3) .eq. 'TIP') resnames(ir)(4:4)='3'
      end do
      if (nconfig .eq. 1) then
        if (ncha .gt. 0) then
          print *,'Atom and/or residue sequence numbers were changed'
        else
          print *,'No sequence number was changed'
        end if
      end if
      return
1000  format(i5)
1001  format(i10)
1002  format(i8)
1006  format(i4)
1100  format(' WARNING: leading digits for ',a,' over ',i7,
     -  '  are dropped')
1101  format(' NOTE:',i5,' atoms with duplicate residue numbers were ',
     -  'inserted',/,6x,'You may want to adjust the residue numbers ',
     -  'to avoid duplicates')
2000  format( 'ERROR: Maximum number of residues (',i5,') is exceeded',
     -  /,' - recompile simulaid with MAXRSD=',i6,' or greater')
      end
