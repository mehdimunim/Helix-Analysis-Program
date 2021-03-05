      subroutine checktorbond(resnami,ixresi,ixresj,atnami,atnamj,
     -  fixbond,ipep,ican,icac,issb)
      character*4 atnami,atnamj
      character*8 resnami,resnam
c     Check for small loops
      fixbond=0
      if (ixresi .eq. ixresj) then
        call leftadjustn(resnami,resnam,8)
        if (resnam(1:3) .eq. 'HIS' .or. resnam(1:3) .eq. 'HSE' .or.
     -      resnam(1:3) .eq. 'HSD' .or. resnam(1:3) .eq. 'HSP' .or.
     -      resnam(1:3) .eq. 'HID') then
          if (atnami(1:3) .eq. 'CD2' .or. atnami(1:3) .eq. 'NE2'
     -    .or. atnami(1:3) .eq. 'CE1' .or. atnami(1:3) .eq. 'ND1')
     -       fixbond=1
          if (atnami(1:3) .eq. 'CG ' .and. atnamj(1:3) .ne. 'CB ')
     -       fixbond=1
        else if (resnam(1:3) .eq. 'PRO') then
          if (atnami(1:3) .eq. 'CD ' .or. atnami(1:3) .eq. 'CB '
     -    .or. atnami(1:3) .eq. 'CG ') fixbond=1
          if ((atnami .eq. 'N   ' .and. atnamj .ne. 'C   ') .or.
     -        (atnami .eq. 'CA  ' .and. atnamj .ne. 'C   '))
     -      fixbond=1
        else if (resnam(1:3) .eq. 'PHE') then
          if (atnami(1:2) .eq. 'CE' .or. atnami(1:2) .eq. 'CD'
     -    .or. atnami(1:2) .eq. 'CZ')  fixbond=1
          if (atnami(1:3) .eq. 'CG ' .and. atnamj(1:3) .ne. 'CB ')
     -       fixbond=1
        else if (resnam(1:3) .eq. 'TYR') then
          if (atnami(1:2) .eq. 'CE' .or. atnami(1:2) .eq. 'CD')
     -      fixbond=1
          if (atnami(1:2) .eq. 'CZ' .and. atnamj(1:2) .ne. 'OH')
     -        fixbond=1
          if (atnami(1:3) .eq. 'CG ' .and. atnamj(1:3) .ne. 'CB ')
     -       fixbond=1
        else if (resnam(1:3) .eq. 'TRP') then
          if (atnami(1:2) .eq. 'CE' .or. atnami(1:2) .eq. 'CD'
     -    .or. atnami(1:2) .eq. 'CZ' .or. atnami(1:2) .eq. 'CH'
     -     .or. atnami(1:2) .eq. 'NE')  fixbond=1
          if (atnami(1:3) .eq. 'CG ' .and. atnamj(1:3) .ne. 'CB ')
     -       fixbond=1
        else if (resnam(1:3) .eq. 'ARG') then
          if (atnami(1:2) .eq. 'NE' .and. atnamj(1:2) .eq. 'CZ' .or.
     -        atnami(1:2) .eq. 'CZ' .and. atnamj(1:2) .eq. 'NH')
     -      fixbond=1
       end if
      end if
      if (fixbond .eq. 0) then
        ipep=0
        ican=0
        icac=0
        issb=0
        if ((atnami .eq. 'C      ' .and. atnamj .eq. 'N      ') .or.
     -    (atnamj .eq. 'C      ' .and. atnami .eq. 'N      ')) ipep=1
        if ((atnami .eq. 'CA     ' .and. atnamj .eq. 'N      ') .or.
     -    (atnamj .eq. 'CA     ' .and. atnami .eq. 'N      ')) ican=1
        if ((atnami .eq. 'CA     ' .and. atnamj .eq. 'C      ') .or.
     -    (atnamj .eq. 'CA     ' .and. atnami .eq. 'C      ')) icac=1
        if (atnami .eq. 'SG     ' .and. atnamj .eq. 'SG     ')  issb=1
       end if
      return
      end
