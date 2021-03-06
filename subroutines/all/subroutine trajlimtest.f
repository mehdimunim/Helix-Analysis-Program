      subroutine trajlimtest(nframe,MAXFRAMES)
      if (nframe .gt. MAXFRAMES) then
        write (6,1004)
        stop
      end if
1004  format(' Calculation stopped since this option is limited to ',
     -  'reading ',i8,' frames',/,' Recompile with larger value of ',
     -  'the parameter MAXFRAMES')
      end
