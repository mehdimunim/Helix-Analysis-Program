      subroutine readint(line,icol1,icol2,intg,ihextyp,istop,ierr)
      character*(*) line
      common /askhex/ iaskhex(4),ishex(4)
      character*1 tab,ctrlM
      common /tab/ tab,ctrlM
      common /logging/ logfile,ipredict
      character*14 inplab(4)
      data inplab /'atom number   ','residue number','residue id    ',
     -  'data          '/
c     ihextyp=1: atom #; ihextyp=2: residue #; inphextyp=4: everything else
      ierr=0
      if (icol1 .gt. icol2) then
        write (6,2001) inplab(ihextyp),icol1,icol2,line(1:icol2)
        stop
      end if
      nchars=0
      do ic=icol1,icol2
        if (line(ic:ic) .ne. ' ' .and. line(ic:ic) .ne. tab)
     -    nchars=nchars+1
      end do
      if (nchars .eq. 0) then
c       Field was blank - return zero
        intg=0
      else if (ishex(ihextyp) .eq. 0) then
        read (line(icol1:icol2),*,ERR=100) intg
      else
c       Interpret string as hexadecimal
        intg=iconvhex(line(icol1:icol2),icol2-icol1+1)
      end if
      return
100   if (iaskhex(ihextyp) .eq. 1) then
c       Check if hexadecimal
        ishex(ihextyp)=ishexadecimal(line(icol1:icol2),icol2-icol1+1)
        if (ishex(ihextyp) .eq. 0) then
          write (6,2000) inplab(ihextyp),icol1,icol2,line
          ierr=1
          if (istop .eq. 1) stop
        else
          print *,'Input ',inplab(ihextyp),' appears to be hexadecimal'
          if (ipredict .eq. 0)
     -      call askyn('Is this correct',15,1,1,ishex(ihextyp),00,0)
          iaskhex(ihextyp)=0
          intg=iconvhex(line(icol1:icol2),icol2-icol1+1)
        end if
      else
        write (6,2000) inplab(ihextyp),icol1,icol2,line
        ierr=1
        if (istop .eq. 1) stop
      end if
      return
2000  format(' ERROR: Invalid integer was read for ',a,/,
     -  8x,'Column range: [',i3,',',i3,':]; the line read:'/,a)
2001  format(' PROGRAM ERROR in readint: icol1 > icol2 (',2i3,')',
     -  ' the line read was:'/,a)
      end
