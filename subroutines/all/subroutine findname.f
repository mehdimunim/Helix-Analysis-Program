      subroutine findname(name,list,ifrst,lenlist,ix,lench)
      character*(*) name,list
      dimension list(lenlist)
      ix=0
c     write (6,*) 'FINDNAME ifrst,lenlist,name=',ifrst,lenlist,name
c      write (77,7711) (list(i),i=1,lenlist)
c7711  format(10a5)
      do i=ifrst,lenlist
        if (name(1:lench) .eq. list(i)(1:lench)) then
          ix=i
          return
        end if
      end do
      return
      end
