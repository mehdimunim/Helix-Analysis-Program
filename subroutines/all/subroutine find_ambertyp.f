      subroutine find_ambertyp(inpt,flag,lflag,form,lform)
      character*(*)flag,form
      character*80 liner
      rewind inpt
      liner(1:1)=' '
      do while (liner(1:lflag) .ne. flag(1:lflag))
        read(inpt,1010,end=100) liner
      end do
      call blankout(liner,1,80)
      read(inpt,1010,end=100) liner
      if (liner(1:7) .eq. '%FORMAT') then
        call lastchar(liner,lc,80)
        form=liner(8:lc)
        lform=lc-7
      end if
      return
100   print *,'Format not found after flag ',flag(1:lflag)
1010  format(a)
      stop
      end
