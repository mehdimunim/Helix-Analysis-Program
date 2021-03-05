      subroutine transform_dist(r_read,r_trans)
      common /transform_dist_dat/ itranstyp,nexp,rijmin,fracexp
      if (itranstyp .eq. 1) then
        r_trans=r_read
      else if (itranstyp .eq. 2) then
        r_trans=r_read**nexp
      else if (itranstyp .eq. 3) then
        r_trans=r_read**fracexp
      else if (itranstyp .eq. 4) then
        r_trans=1.0-r_read
      else if (itranstyp .eq. 5) then
        r_trans=(1.0-r_read)**nexp
      else if (itranstyp .eq. 6) then
        r_trans=(1.0-r_read)**fracexp
      else if (itranstyp .eq. 7) then
        r_trans=abs(r_read)**nexp
      else if (itranstyp .eq. 8) then
        r_trans=(1.0-abs(r_read))**nexp
      else if (itranstyp .eq. 9) then
        r_trans=(r_read-rijmin)**nexp
      end if
      return
      end
