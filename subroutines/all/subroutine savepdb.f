      subroutine savepdb(iu,label,llabel,c,n,mod)
      character*(*) label
      dimension c(3,n)
      parameter (MAXREC=200000)
      character*132 line(MAXREC)
      dimension index(MAXREC)
      common /line_crd/ line,index
      if (mod .le. 1) then
        open(unit=iu,status='new',file=label(1:llabel),iostat=iopen,
     -    form='formatted')
        if (iopen .gt. 0)
     -    open(unit=iu,status='old',file=label(1:llabel),iostat=iopen,
     -        form='formatted')
        if (iopen .eq. 0) then
          print *,'File ',label(1:llabel),' opened'
        else
          return
        end if
        write (iu,1002) label(1:llabel)
      end if
      if (mod .gt. 0) write (iu,1003) mod
      do ia=1,n
        write (line(index(ia))(31:54),1001) (c(k,ia),k=1,3)
        call lastchar(line(index(ia)),lc,132)
        write (iu,1000) line(index(ia))(1:lc)
      end do
      if (mod .gt. 0) then
        write (iu,1004)
      else
        close(iu)
      end if
      return
1000  format(a)
1001  format(3f8.2)
1002  format('REMARK ',a)
1003  format('MODEL ',i4)
1004  format('ENDMDL')
      end
