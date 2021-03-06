      subroutine readcharges(nread,nslt,n,charge,iatnum,isv,
     -  icharges,nerr)
      dimension charge(nread),iatnum(nread)
      character*1 ans
      character*8 atomnam
      character*80 line
      character*200 filename
      inpt=51
      call quiz(ans,ians,'n',' ',0,'charge input',12,0,5,6,0)
      if (ans .eq. 'n') return
      nread=nslt
      if (n .gt. nslt) then
        call askyn('Do you have charges for the solvents too',40,1,1,
     -    isv,000,0)
        if (isv .eq. 1) nread=n
      end if
      namlen=0
      if (ans .eq. 'a') then
        call openfile(inpt,0,'Amber prtmtop',13,'old',filename,namlen,
     -    notfound,3,1,1,0,0)
        do while (line(1:12) .ne. '%FLAG CHARGE')
          read (inpt,1000,end=199) line
        end do
        read (inpt,1000,end=199) line
        read (inpt,1001,end=198,err=198) (charge(i),i=1,nread)
        do i=1,nread
          charge(i)=charge(i)/18.2223
c         write (77,*) i,' q=',charge(i)
        end do
      else if (ans .eq. 'c') then
        call openfile(inpt,0,'Charmm PSF',10,'old',filename,namlen,
     -    notfound,3,1,1,0,0)
        call find_n_psf(inpt,line,nerr,nread,natspsf,'!NATOM',6)
        if (nerr .gt. 0) go to 999
        do ia=1,nread
          call blankout(line,1,80)
          read (inpt,1000,end=999) line
          ic=1
          do irec=1,7
            call nextstring(line,ic,ic1,ic2,80)
            if (irec .eq. 5) then
              atomnam='        '
              atomnam(1:ic2-ic1+1)=line(ic1:ic2)
            end if
          end do
          charge(ia)=0.0
          iok=0
          read (line(ic1:ic2),*,err=603) charge(ia)
          iok=1
603       if (iok .eq. 0) then
            write (6,2000) 'PSF',line(1:79)
            nerr=nerr+1
           end if
        end do
      else if (ans .eq. 'd') then
        call openfile(inpt,0,'Autodock .pdbqt',15,'old',filename,namlen,
     -    notfound,3,1,1,0,0)
        ia=0
        do while (ia .lt. nread)
          read (inpt,1000,end=399) line
          if (line(1:4) .eq. 'ATOM' .or. line(1:6) .eq. 'HETATM') then
            ia=ia+1
            charge(ia)=0.0
            iok=0
            read (line(71:76),*,err=604) charge(ia)
            iok=1
604         if (iok .eq. 0) then
              write (6,2000) '.pdbqt',line(1:79)
              nerr=nerr+1
            end if
          end if
        end do
399     if (ia .lt. nread) then
          nerr=1
          write (6,2001) ia,nread
        end if
      end if
      nwarn=0
      qsum=0.d0
      nqz=0
      do i=1,nread
        qsum=qsum+charge(i)
        if (charge(i) .eq. 0.0) nqz=nqz+1
        if (charge(i) .lt. -1.0001 .or. charge(i) .gt. 1.0001) then
          nwarn=nwarn+1
          write (6,2004) i,charge(i)
        else if (iatnum (i) .gt. 0) then
          if (iatnum(i) .eq. 1 .and. charge(i) .lt. -0.05 .or.
     -        iatnum(i) .eq. 8 .and. charge(i) .gt. 0.0) then
            write (6,2002) i,iatnum(i),charge(i)
            nwarn=nwarn+1
          end if
        end if
      end do
      if (nqz .eq. nread) then
        print *,'All charges read are zero'
        nerr=nerr+1
      else
        write (6,2003) qsum,nqz
      end if
      if (nerr .gt. 0) then
        print *,'There were errors in reading the charges'
        call askstop(1)
        print *,'No charges will be used'
      end if
      go to 999
199   print *,'Did not find charge flag in file ',filename(1:namlen)
      nerr=1
      go to 999
198   print *,'Run out of data'
      nerr=1
999   close (inpt)
      if (nerr .eq. 0) icharges=2
      return
1000  format(a)
1001  format(5e16.8)
2000  format(' ERROR: invalid charge in ',a,' record:',/,1x,a)
2001  format(' ERROR: only ',i7,' atoms were found (instead of ',i7,')')
2002  format(' WARNING: Atom',i6,' atomic #',i2,' has unlikely sign:',
     -  f5.2)
2003  format(' Total charge =',f8.3,/,' Number of zero charges=',i6)
2004  format(' WARNING: charge for atom',i6,': ',f8.4)
      end
