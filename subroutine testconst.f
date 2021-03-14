      subroutine testconst(izero,ione,itwo,fzero,fone,ftwo,iout,nfail,
     -  initfail,lab)
c#    MMC routine 105 lstmod: 01/18/07
      character*4 lab
c*****Check the integrity of stored constants
      if (initfail .gt. 0) nfail=0
      if (itwo+itwo .ne. itwo*itwo .or. itwo+itwo .eq. itwo) then
        nfail=nfail+1
        write (iout,1000) lab,'integer two'
      end if
      if (ftwo+ftwo .ne. ftwo*ftwo .or. ftwo+ftwo .eq. ftwo) then
        nfail=nfail+1
        write (iout,1000)lab, 'float two'
      end if
      if (ione+ione .eq. ione .or. ione*ione .ne. ione) then
        nfail=nfail+1
        write (iout,1000) lab,'integer one'
      end if
      if (fone+fone .eq. ione .or. fone*fone .ne. fone) then
        nfail=nfail+1
        write (iout,1000) lab,'float one'
      end if
      if (izero+izero .ne. izero .or. izero*izero .ne. izero) then
        nfail=nfail+1
        write (iout,1000) lab,'integer zero'
      end if
      if (fzero+fzero .ne. fzero .or. fzero*fzero .ne. fzero) then
        nfail=nfail+1
        write (iout,1000) lab,'float zero'
      end if
      return
1000  format(' ***** PROGRAM ERROR at ',a,': Internal constant ',a,
     -  ' is corrupted')
      end
