      subroutine changeext(oldname,newname,leno,lenn,newext,lenext,
     -  insert,iask)
      character*(*) oldname,newname,newext
      character*80 question,oldext
      character*200 tempname
      tempname=oldname
      lenn=leno
      do while (lenn .gt. 0 .and. tempname(lenn:lenn) .ne. '.')
        lenn=lenn-1
      end do
      if (lenn .eq. 0) then
c       Input file had no '.' in it - just add extension
        lenn=leno+1
        tempname(lenn:lenn)='.'
      end if
      if (insert .eq. 0) then
c       Replace the extension
        tempname(lenn+1:lenn+lenext+1)=newext
        lenn=lenn+lenext+1
        if (iask .eq. 0) then
          newname=tempname
        else
          question='Output file name: '//tempname(1:lenn)//' - is it OK'
          lq=29+lenn
          call askyn(question,lq,1,1,ians,00,0)
          if (ians .eq. 1) then
            newname=tempname
          else
            lenn=0
          end if
        end if
      else
c       Insert extension between file name root and old extension
        oldext=tempname(lenn:leno)
        newname=tempname(1:lenn-1)//'_'//newext//oldext
        lenn=leno+lenext+1
      end if
      return
      end
