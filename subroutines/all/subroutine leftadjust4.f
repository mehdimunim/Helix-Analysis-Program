      subroutine leftadjust4(in,out)
c*****Left-adjust a string of four characters
      character*(*) in
      character*4 out,outt
      outt='    '
      do i=1,4
        if (in(i:i) .ne. ' ') then
          do j=i,4
            outt(j-i+1:j-i+1)=in(j:j)
          end do
          out=outt
          return
        end if
      end do
      out=in(1:4)
      return
      end
