      subroutine unpackconf(line,nlines,ncol,inpfile,namleni,nconfig,
     -  molname,molnamelen,numsel,iconfsel,numconsec,maxconfsel,
     -  nextconfsel,incr_fileno,ioverallconf,ionefile,
     -  natsfirst,nats,iuunp,maxrec)
c     Unpack the next conformation from the stack
      dimension iconfsel(maxconfsel)
      character*(*) molname
      character*132 line(maxrec)
      character*200 inpfile,outfile
      common /logging/ logfile,ipredict
c     print *,'UNPACKCONF nextconfsel,nconfig,ionefile=',
c    -                    nextconfsel,nconfig,ionefile
      if (numsel .gt. maxconfsel) then
        print *,'PROGRAM ERROR: illegal numsel (',numsel,');',
     -    ' maxconfsel=',maxconfsel
        stop
      end if
      if (numsel .gt. 0) then
        if (nconfig .eq. iconfsel(nextconfsel)) then
c         Selection found - increment nextconfsel
          nextconfsel=nextconfsel+1
        else
c         Configuration not selected - skip
          return
        end if
      end if
c     print *,'Unpackconf nconfig,numsel,nextconfsel=',
c    -  nconfig,numsel,nextconfsel
      if (nconfig .eq. 1) natsfirst=nats
      if (nextconfsel .le. 10) then
        if (nats .ne. natsfirst) write (6,1012) nconfig,nats,natsfirst
      end if
      if (ionefile .eq. 0) then
        if (molnamelen .gt. 0) then
c         Currently only DOCK PDB files have molecule names read
          nl1=molnamelen
          outfile(1:nl1)=molname(1:nl1)
          outfile(nl1+1:nl1+4)='.pdb'
          nl1=nl1+4
        else
          ifnumw=nconfig+incr_fileno
          if (numconsec .gt. 0) ifnumw=nextconfsel-1+incr_fileno
          call filenamenum(inpfile,namleni,outfile,nl1,ifnumw,+2)
        end if
        call openfile(iuunp,0,'output',6,'new',outfile,nl1,notfnd,2,1,1,
     -    0,ioverallconf)
        if (notfnd .eq. 1) write (6,1011) outfile(1:nl1)
        if (notfnd .eq. 1 .or. 
     -      (ipredict .eq. 1 .and. ioverallconf .eq. 0)) then
          call askyn('Do you want to overwrite all existing files',43,
     -      1,-1,ioverallconf,0,0)
          if (ioverallconf .eq. 0) stop
          call openfile(iuunp,0,'output',6,'new',outfile,nl1,notfnd,2,1,
     -      1,0,ioverallconf)
        end if
        write (6,1010) outfile(1:nl1)
      else if (nextconfsel .eq. 2) then
c       Open overall file (input filename with_sel added to it)
        outfile=inpfile
        ic=namleni
        do while (ic .gt. 1 .and. inpfile(ic:ic) .ne. '.')
          ic=ic-1
        end do
        if (ic .gt. 1) then
          outfile(ic:ic+3)='_sel'
          outfile(ic+4:namleni+4)=inpfile(ic:namleni)
        else
          outfile(namleni+1:namleni+4)='_sel'
        end if
        nl1=namleni+4
        call openfile(iuunp,0,'output',6,'new',outfile,nl1,notfnd,0,1,
     -    1,0,ioverallconf)
        write (6,1010) outfile(1:nl1)
      end if
      lf=1
      if (line(1)(1:3) .eq. 'TER' .or. line(1)(1:3) .eq. 'END') lf=2
      do i=lf,nlines
        call writeline(iuunp,line(i),1,ncol,0)
      end do
      if (ionefile .eq. 0) close (iuunp)
      return
1010  format(' Writing file: ',a)
1011  format(' Problem opening file: ',a)
1012  format(' NOTE: # of atoms in configuration',i7,'=',i6,
     -  'in the first conf:',i6)
      end
