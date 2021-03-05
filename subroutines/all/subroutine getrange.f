      subroutine getrange(ifst,ifstdef,ilst,ilstdef,incr,incrask,
     -  question,lquestion,maxval,ihelp)
      character*(*) question
      character*80 line
c     Input a range
c     print *,'GETRANGE ifst,ifstdef,ilst,ilstdef=',
c    - ifst,ifstdef,ilst,ilstdef
c     print *,'GETRANGE maxval=',maxval
      line(7:6+lquestion)=question(1:lquestion)
100   line(1:6)='First '
      call getint(line,6+lquestion,ifstdef,1,maxval,ifst,ihelp)
      line(1:6)='Last  '
      call getint(line,6+lquestion,ilstdef,1,maxval,ilst,ihelp)
      if (ilst .lt. ifst) then
        print *,'Invalid range'
        go to 100
      end if
      if (incrask .gt. 0) then
101     call getint('Increment',9,1,1,ilst-ifst+1,incr,47)
        if (incr .lt. 1) then
          print *,'Invalid increment'
          go to 101
        end if
      end if
      return
      end
