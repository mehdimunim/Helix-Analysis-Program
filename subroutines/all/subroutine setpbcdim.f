      subroutine setpbcdim(ioppbc,ixyzhex,ixyzexcld,ixyzincld,xyz)
      dimension ixyzhex(3)
      character*1 xyz(3),ans
      ixyzexcld=0
      ixyzincld=0
      if (ioppbc .ne. 3 .and. ioppbc .ne. 5 .and.
     -    ioppbc .ne. 6 .and. ioppbc .ne. 7 .and.
     -    ioppbc .ne. 9 .and. ioppbc .ne. -1) then
        call getint('Number of PBC dimensions',24,3,1,3,npbcdim,110)
        if (npbcdim .eq. 2) then
          if (ioppbc .eq. 4 .or. ioppbc .eq. 5) then
            ixyzexcld=ixyzhex(1)
            print *,'Hexagonal prism - exluded axis:',xyz(ixyzexcld)
          else
            call quiz(ans,ixyzexcld,' ',
     -        'to exclude from the PBC calculation',-35,'axis',4,
     -        0,5,6,0)
          end if
        else if (npbcdim .eq. 1) then
          if (ioppbc .eq. 4 .or. ioppbc .eq. 5) then
            ixyzincld=ixyzhex(1)
            print *,'Hexagonal prism - included axis:',xyz(ixyzincld)
          else
            call quiz(ans,ixyzincld,' ',
     -        'to include in the PBC calculation',-33,'axis',4,
     -         0,5,6,0)
          end if
        end if
      end if
      return
      end
