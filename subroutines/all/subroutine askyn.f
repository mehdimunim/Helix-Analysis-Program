      subroutine askyn(q,lenqin,iyn,idefans,ians,ihelp,itip)
      character*(*) q
      character*132 pline
      character*1 ans
      character*5 defans
      common /graphics/ npixx,npixy,maxpixx,maxpixy,idwmain,idwplot,
     -  wx,wy,wz,wxdr
      common /logging/ logfile,ipredict
c     idefans=-1: default no; idefans=+1: default yes
c     iyn=1: yes -> ians=1, no -> ians=0
c     iyn=0: yes -> ians=0, no -> ians=1
c     print *,'ASKYN ihelp,itip=',ihelp,itip
      ians=0
      if (idefans .eq. -1) then
        defans=' [n] '
        lendef=5
      else if (idefans .eq. +1) then
        defans=' [y] '
        lendef=5
      else
        defans=' '
        lendef=1
      end if
      lenq=lenqin
      pline(1:1)=' '
      pline(2:lenq+1)=q(1:lenq)
      pline(lenq+2:lenq+6)=' (y/n'
      lenq=lenq+6
      if (ihelp .gt. 0) then
        pline(lenq+1:lenq+2)='/?'
        lenq=lenq+2
      end if
      if (itip .gt. 0) then
        pline(lenq+1:lenq+2)='/+'
        lenq=lenq+2
      end if
      lenq=lenq+1
      pline(lenq:lenq)=')'
      pline(lenq+1:lenq+lendef)=defans(1:lendef)
100   call writeline(6,pline,1,lenq+lendef,1)
      ans= ' '
      read (5,1001,end=99,err=99) ans
      if (ans .eq. '?') then
        if (ihelp .eq. 0) then
          print *,'Sorry, no help on this'
          go to 100
        else
          call explanation(ihelp,0)
          go to 100
        end if
      end if
      if (ans .eq. '+') then
        if (itip .eq. 0) then
          print *,'Sorry, no tip on this'
          go to 100
        else
          call explanation(0,itip)
          go to 100
        end if
      end if
99    if (ans .ne. 'n' .and. ans .ne. 'N' .and. ans .ne. 'y' .and.
     -  ans .ne. 'Y') then
        call lastchar(ans,ilc,1)
        if (ilc .eq. 1) then
          print *,'Invalid answer - please, answer y or n'
          go to 100
        else if (idefans .ne. 1 .and. idefans .ne. -1) then
          print *,'Pls answer y or n'
          go to 100
        else
          if (idefans .eq. -1) then
            ans='n'
          else if (idefans .eq. 1) then
            ans='y'
          end if
        end if
      end if
      if (logfile .gt. 0) write (logfile,1001) ans
      if (ans .eq. 'y' .or. ans .eq. 'Y') ians=1
      if (iyn .eq. 0) ians=1-ians
      return
1001  format(a1)
      end
