      subroutine addatsym(line,ian,ifc)
      character*(*) line
      character*2 iatnm2
      common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
     -  mmatno(64),iatnm2(99)
      character*2 chemsym
      character*4 name
c     print *,'ADDATSYM ifc=',ifc
c     Add atomic symbol to col 77-78
      if (ian .lt. 1) then
        name=line(13:16)
        iatno=ianum(name,1,4)
      else
        iatno=ian
      end if
      chemsym=iatnm2(iatno)
      if (chemsym(2:2) .eq. ' ') then
         chemsym(2:2)=chemsym(1:1)
         chemsym(1:1)=' '
      end if
      if (ifc .gt. 76 .and. line(77:78) .ne. chemsym) then
        write (6,1000) line(77:78),chemsym,line(1:ifc)
      else
        call blankout(line,ifc+1,80)
      end if
      line(77:78)=chemsym
      if (ifc .lt. 77) ifc=78
      return
1000  format(' WARNING: columns 77-78 is not blank:',a,' - changed to ',
     -  a,' in line',/,a)
      end
