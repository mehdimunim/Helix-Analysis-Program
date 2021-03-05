      subroutine find_n_psf(inpt,line,ierr,n,n_psf,label,llabel)
      character*80 line
      character*(*) label
      ierr=0
      n_psf=0
      do while (n_psf .eq. 0)
        call blankout(line,1,80)
        read (inpt,1010,end=602) line
        ic=1
        call nextstring(line,ic,ic1,ic2,80)
        call nextstring(line,ic,ic21,ic22,80)
        if (ic22 .gt. ic21) then
          if (line(ic21:ic21+llabel-1) .eq. label(1:llabel))
     -      read(line(ic1:ic2),*,err=602) n_psf
        end if
      end do
602   if (n_psf .eq. 0) then
        print *,'ERROR: no ',label(1:llabel),' line was found'
        ierr=1
      else if (n_psf .lt. n) then
        write (6,2000) n,n_psf
        ierr=1
      else
        print *,'Number of ',label(2:llabel),' found in the PSF file=',
     -    n_psf
      end if
      return
1010  format(a)
2000  format(' ERROR: the number of atoms to read (',i8,') exceeds ',
     -  'the number of atoms in the PSF file (',i8,')')
      end
