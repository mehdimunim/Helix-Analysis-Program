      subroutine checkhbclose(c,n,i1,i2,label,llabel,ihbdel)
      dimension c(3,n)
      character*(*) label
      dij2=dist2(c(1,i1),c(1,i2))
      if (dij2 .lt. 0.5) then
        if (dij2 .lt. 0.001) then
          write (6,1000) 'ERROR',i1,i2,label(1:llabel),sqrt(dij2)
          ihbdel=1
        else
          write (6,1000) 'WARNING',i1,i2,label(1:llabel),sqrt(dij2),
     -      ' - hydrogen bond is deleted'
        end if
      end if
      return
1000  format(1x,a,': atoms ',i6,' and ',i6,' are listed as ',a,/,
     -  ' but are only',f8.5,' A apart',a)
      end
