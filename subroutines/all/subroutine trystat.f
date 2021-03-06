      subroutine trystat(opt,try,ntry,label,llabel)
      dimension try(200)
      character*(*) label
      nt=min0(200,ntry)
      nopt=0
      do it=1,nt
        if (abs(try(it)-opt)/abs(abs(try(it))+abs(opt)) .lt. 0.00001)
     -    nopt=nopt+1
      end do
      write (6,1000) label(1:llabel),nopt
      if (nopt .eq. 0) print *,'WARNING: best orientation was ',
     -  'optimzed from the input conformation'
      if (nopt .eq. 1) print *,
     -  'WARNING: more random trys are likely to reach better optimum'
      return
1000  format(' The number of times the optimal ',a,' was reached=',i3)
      end
