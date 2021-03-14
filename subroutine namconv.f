      subroutine namconv(nrescol,resnamin,atnamin,resnamout,atnamout,
     -  nrch,nach,line,idcol,nrecdel)
c     Perform an atom and residue name conversion
      character* 132 line
      character*1 ans1
      character*4 atnamin,atnamout,an0,an1,resnamin,resnamout,rn0
      character*8 namnam
      common /convspec/ incon,ioutcon,ideoxy,namnam(10)
      character*3 ires0
      character*4 anamcontab,rnamcontab
      common /convdat/ nanamcon,nrnamcon,nres0,ifst(100),ilst(100),
     -  ncon,nconcol,anamcontab(1000,13),rnamcontab(13,100),ires0(100),
     -  icongro
      data inamconv /0/
      incol1=incon+1
      incol2=incon+1
      if (incon .eq. ncon) incol2=nconcol
      an0=atnamin
      rn0=resnamin
c     write (77,*) 'Aresnamin=',resnamin,' atnamin=',atnamin
      atnamout=atnamin
      call leftadjust4(resnamin,resnamin)
      call leftadjust4(atnamin,an1)
      inrescol=nrescol
      if (incon .eq. ncon) then
        call regularpdb(an1,atnamin,1)
        inrescol=3
      else
        if (idigit(an1(1:1),1) .eq. 1 .and. ioutcon .ne. ncon)
     -    call regularpdb(an1,atnamin,-1)
      end if
c     write (77,*) ' an1=',an1,' atnamin=',atnamin
c     print *, 'an1,atnamin=',an1,'|',atnamin
      atnamout=atnamin
      an1=atnamin
      resnamout=resnamin
c     Find first the residue name in the resname table
      iresconv=0
      do i=1,nrnamcon
        do j=incol1,incol2
          if (resnamin(1:inrescol) .eq.
     -        rnamcontab(j,i)(1:inrescol)) then
            iresconv=i
            resnamout=rnamcontab(ioutcon+1,i)
c           write (77,*) 'Match i,j=',i,j,' rno=',resnamout
            go to 100
          end if
        end do
      end do
      go to 910
c     Find generic residue name
100   iresconv0=0
      do i=1,nres0
        if (ires0(i) .eq. rnamcontab(1,iresconv)(1:3)) then
          iresconv0=i
c         write (77,*) 'Generic resname=',ires0(i)
          go to 110
        end if
      end do
      write (6,2001) rnamcontab(1,iresconv)
c     Find the atom conversion - first check residue-specific conversions
110   if (iresconv0 .gt. 0) then
        do i=ifst(iresconv0),ilst(iresconv0)
          do j=incol1,incol2
c           write (77,*)
c    -        'i=',i,' atnamin=',atnamin,' atab=',anamcontab(i,j)
            if (atnamin .eq. anamcontab(i,j)) then
              if (anamcontab(i,ioutcon+1) .eq. '****') then
                line(idcol:idcol)='*'
                nrecdel=nrecdel+1
              else if (anamcontab(i,ioutcon+1) .ne. '    ') then
                atnamout=anamcontab(i,ioutcon+1)
              end if
c             write (77,7712) atnamin,resnamin,atnamout
c7712          format(' FOUND: an=',a4,' rn=',a4,' ao=',a4)
              inamconv=i
              go to 900
            end if
          end do
        end do
      end if
c     Not found among the specifics  - check the residue independent list
      do i=ifst(1),ilst(1)
        do j=incol1,incol2
          if (atnamin .eq. anamcontab(i,j)) then
            if (anamcontab(i,ioutcon+1) .eq. '****') then
              line(idcol:idcol)='*'
              nrecdel=nrecdel+1
            else if (anamcontab(i,ioutcon+1) .ne. '    ') then
              atnamout=anamcontab(i,ioutcon+1)
            end if
            inamconv=i
c           write (77,*) 'Generic found atnamout=',atnamout,' i,j=',i,j
            go to 900
          end if
        end do
      end do
      inamconv=0
c     Take care of oxy/deoxy nucleic acids
900   if (nrescol .ge. 4 .and. resnamout(4:4) .ne. ' ') then
        if (resnamout(1:1) .eq. 'D' .or. resnamout(1:1) .eq. 'R') then
          if (ideoxy .eq. -1) then
            if (resnamout(2:4) .eq. 'ADE' .or. resnamout(2:4) .eq.
     -        'GUA' .or. resnamout(2:4) .eq. 'CYT') then
c             Ask if oxy or deoxy NA
              call getname(ans1,len,
     -          'Nucleic acid type (D for Deoxy, R Oxy)',38,1,'',0,0,0,
     -          0)
              ideoxy=1
              if (ans1 .eq. 'r' .or. ans1 .eq. 'R') ideoxy=0
            else if (resnamout(2:4) .eq. 'THY') then
              ideoxy=1
            else if (resnamout(2:4) .eq. 'URA') then
              ideoxy=0
            end if
          else if (ideoxy .eq. 1) then
            resnamout(1:1)='D'
          else if (ideoxy .eq. 0) then
            resnamout(1:1)='R'
          end if
        end if
      end if
      if (ioutcon .eq. icongro) then
c       Gromacs output - check for residue name change
        ifound=0
        if (ires0(iresconv0) .eq. 'ASP') then
          call findname('HD2 ',anamcontab(1,ioutcon+1),ifst(iresconv0),
     -      ilst(iresconv0),ifound,3)
        else if (ires0(iresconv0) .eq. 'CYS') then
          call findname('HG  ',anamcontab(1,ioutcon+1),ifst(iresconv0),
     -      ilst(iresconv0),ifound,2)
        else if (ires0(iresconv0) .eq. 'GLU') then
          call findname('HE2 ',anamcontab(1,ioutcon+1),ifst(iresconv0),
     -      ilst(iresconv0),ifound,3)
        else if (ires0(iresconv0) .eq. 'LYS') then
          call findname('HZ3 ',anamcontab(1,ioutcon+1),ifst(iresconv0),
     -      ilst(iresconv0),ifound,3)
        end if
        if (ifound .gt. 0) resnamout(4:4)='H'
      end if
c     Check for changes
910   if (inamconv*iresconv .eq. 0) then
        if (incon .eq. ncon) then
c         Undo the regularization
          call leftadjust4(an0,an1)
          if (idigit(an1(1:1),1) .eq. 1)
     -      call regularpdb(an1,atnamout,-1)
        else
          atnamout=an1
        end if
      end if
      if (rn0 .ne. resnamout) nrch=nrch+1
      if (an0 .ne. atnamout) nach=nach+1
c      write (77,7723) rn0,an0,resnamout,atnamout
c7723  format(' rn0,an0=',a,1x,a,' ,resnamout,atnamout=',a,1x,a)
      return
2001  format(' ERROR: Generic residue id ',a4,' is not found - ',
     -  'conversion is skipped',/,8x,'Check the conversion table file')
      end
