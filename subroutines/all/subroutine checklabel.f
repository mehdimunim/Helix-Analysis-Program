      subroutine checklabel(line,label,llabel,inpt)
      character*1000 line
      character*10 label,labelread
      call blankout(line,1,80)
      read (inpt,1000) line
      call lastchar(line,ilc,80)
c     write (77,*) 'CHECKLABEL',line(1:ilc)
      llabelread=ilc-4
      labelread(1:llabelread)=line(5:ilc)
      if (labelread(1:llabelread) .ne. label(1:llabel)) then
        print *,'ERROR: incorrect data label:',labelread(1:llabelread),
     -    ' (instead of ',label(1:llabel),')'
        stop
      end if
      return
1000  format(a)
      end
