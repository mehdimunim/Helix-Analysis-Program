      subroutine progress_rep(nframe,nframe2d,nframesign)
      nfr=max0(nframe,nframe2d)
      if (mod(nfr,nframesign) .eq. 0) then
        iperc=nfr/nframesign
        if (iperc .lt. 10) write (6,8099) iperc
      end if
8099  format(' Trajectory scan',i3,'0% done')
      end
