      subroutine length_test_check(lentest_ok,iout)
      dimension icntrl(21)
      lentest_ok=1
      call indexit(icntrl,1,20,0)
      open(unit=10,status="new",file='intel_test',form="unformatted",
     -  iostat=iopenok)
      if (iopenok .gt. 0) open(unit=10,status="old",file='intel_test',
     -   form="unformatted",err=100)
      write (10) (icntrl(i),i=1,20)
      write (10) 100
      close (10)
      open(unit = 10 , status = "OLD" , file = 'intel_test',
     -   form="UNFORMATTED" , ERR = 100)
      icntrl(21)=-999
      read(10,end=200,err=200) icntrl
200   if (icntrl(21) .ne. -999) then
        write (iout,*) 'NOTE: record length checks will be skipped'
        lentest_ok=0
      end if
      close (10,status='delete')
      return
100   print *,'PROGRAM ERROR: could not open file intel_test'
      return
      end
