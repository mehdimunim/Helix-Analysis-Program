      subroutine findmenu(label,llabel,ix)
      character*(*) label
      character*60 promptlist,prompttype
      common /quizinfo/ nqst(600),lqst(600),iqfst(600),iqlst(600),
     -  maxq,init,nprompttype,promptlist(600),prompttype(600)
      ix=1
      do while (promptlist(ix)(1:5) .ne. '**** ' .or.
     -          promptlist(ix)(6:5+llabel) .ne. label(1:llabel))
        ix=ix+1
        if (ix .eq. 500) then
          print *,'PROGRAM ERROR: menu ',label(1:llabel),' is not found'
          stop
        end if
      end do
      return
      end
