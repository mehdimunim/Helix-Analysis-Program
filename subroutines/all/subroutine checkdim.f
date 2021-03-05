      subroutine checkdim(n,max,param,lparam,name,lname,iout)
      character*(*) param,name
      if (n .gt. max) then
        write (6,2000) name(1:lname),n,max,param(1:lparam)
        if (iout .gt. 0)
     -    write (iout,2000) name(1:lname),n,max,param(1:lparam)
        stop
      end if
      return
2000  format(' ERROR: ',a,' (',i8,') exceeds limit (',i8,')',/,
     -  ' Increase the parameter ',a,' and recompile')
      end
