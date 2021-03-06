      subroutine updatesolvents(iaskatnum,iaskname,nmax,n,naslv,iasv,
     -  namesv,resnamslv,qsv,pflsv,iwat,iatnum,nconf,iotyp,maxrec)
      character*4 pflsv(100)
      character*8 namesv,resnamslv
      dimension iasv(100),namesv(100),qsv(100),iatnum(maxrec)
      character*5 crdext
      character*80 question
      common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
     -  iommod,iommc,iommc4,iogro,iomol2,iomae,iocif,ioxxx,ioins,
     -  ionxyz,iosxyz,iosxyzrq,iograsp,iofull,lext(19),crdext(19)
c     print *,'UPDATESOLVENTS nmax,n,naslv=',nmax,n,naslv
      if (iaskatnum .eq. 0) then
        call askyn('Is the solvent water (O,H,H)',28,1,1,iwat,0,0)
        if (iwat .eq. 1) then
          naslv=3
          iasv(1)=8
          iasv(2)=1
          iasv(3)=1
          namesv(1)='OH2  '
          namesv(2)='H1   '
          namesv(3)='H2   '
          call getreal('Oxygen charge',13,-0.834,qox,0,0)
          qsv(1)=qox
          qsv(2)=-qox/2.0
          qsv(3)=-qox/2.0
          if (iotyp .eq. ioins .or. iotyp .eq. iommc) then
            if (iotyp .eq. iommc)
     -        print *,"TIP3P oxygen is 'OT' and hydrogen is 'HT'"
            if (iotyp .eq. ioins)
     -        print *,"CVFF oxygen is 'o*' and hydrogen is 'h*'"
            pflsv(1)='    '
            question='Atomtype label for water oxygen'
            call getname(pflsv(1),len,question,31,4,'',0,0,0,0)
            pflsv(2)='    '
            question(26:32)='hydrogen'
            call getname(pflsv(2),len,question,32,4,'',0,0,0,0)
            pflsv(3)=pflsv(2)
          end if
        else
          call getint('Number of atoms in a solvent molecule',37,
     -      999999,1,0,naslv,0)
          do i=1,naslv
            question='Atomic number for solvent atom   '
            write (question(31:33),1000) i
            call getint(question,33,999999,1,99,iasv(i),00)
            if (iaskname .eq. 1) then
              question='Atom name for solvent atom   '
              write (question(27:29),1000) i
              call getname(namesv(i),len,question,29,5,'',0,0,0,0)
              iasv(i)=ianum(namesv(i),0,len)
              call getreal('Charge',6,999999.0,qsv(i),0,0)
              if (iotyp .eq. ioins .or. iotyp .eq. iommc) then
              question='Atom type label for solvent atom   '
              write (question(33:35),1000) i
              call getname(pflsv(i),len,question,35,4,'',0,0,0,0)
              end if
            end if
          end do
        end if
        if (iaskname .eq. 1 .and. resnamslv .eq. '     ')
     -    call getname(resnamslv,len,'Residue name of the solvent',27,
     -      5,'',0,0,0,0)
        iaskatnum=1
      end if
      if (mod(n-nmax,naslv) .ne. 0) then
        if (nconf .lt. 10) then
          write (6,1003) naslv,n,nmax
          call askstop(0)
        end if
      end if
      ia=nmax
      nmax=n
      do while (ia .lt. n)
        do ias=1,naslv
          iatnum(ia+ias)=iasv(ias)
        end do
        ia=ia+naslv
      end do
      return
1000  format(i3)
1003  format(' WARNING: extra atoms dont form full solvents',/,
     -  ' Number of atoms/solvent=',i5,/,
     -  ' Total number of atoms=',i10,/,
     -  ' Number of already existing atoms=',i10)
      end
