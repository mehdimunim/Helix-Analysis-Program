      subroutine readreal(line,icol1,icol2,r)
      character*(*) line
      if (icol2-icol1 .eq. 2) then
        read (line(icol1:icol2),1003,ERR=100) r
      else if (icol2-icol1 .eq. 3) then
        read (line(icol1:icol2),1004,ERR=100) r
      else if (icol2-icol1 .eq. 4) then
        read (line(icol1:icol2),1005,ERR=100) r
      else if (icol2-icol1 .eq. 5) then
        read (line(icol1:icol2),1006,ERR=100) r
      else if (icol2-icol1 .eq. 6) then
        read (line(icol1:icol2),1007,ERR=100) r
      else if (icol2-icol1 .eq. 7) then
        read (line(icol1:icol2),1008,ERR=100) r
      else if (icol2-icol1 .eq. 8) then
        read (line(icol1:icol2),1009,ERR=100) r
      else if (icol2-icol1 .eq. 9) then
        read (line(icol1:icol2),1010,ERR=100) r
      else if (icol2-icol1 .eq. 10) then
        read (line(icol1:icol2),1011,ERR=100) r
      else if (icol2-icol1 .eq. 11) then
        read (line(icol1:icol2),1012,ERR=100) r
      else if (icol2-icol1 .eq. 12) then
        read (line(icol1:icol2),1013,ERR=100) r
      else if (icol2-icol1 .eq. 13) then
        read (line(icol1:icol2),1014,ERR=100) r
      else if (icol2-icol1 .eq. 14) then
        read (line(icol1:icol2),1015,ERR=100) r
      else if (icol2-icol1 .eq. 15) then
        read (line(icol1:icol2),1016,ERR=100) r
      else if (icol2-icol1 .eq. 16) then
        read (line(icol1:icol2),1017,ERR=100) r
      else if (icol2-icol1 .eq. 17) then
        read (line(icol1:icol2),1018,ERR=100) r
      else if (icol2-icol1 .eq. 18) then
        read (line(icol1:icol2),1019,ERR=100) r
      else if (icol2-icol1 .eq. 19) then
        read (line(icol1:icol2),1020,ERR=100) r
      else if (icol2-icol1 .gt. 20) then
        print *,'EXTEND subroutine readreal to larger reals'
        print *,'r=',r,' icol1,2=',icol1,icol2
        stop
      else
        print *,'PROGRAM ERROR: invalid column range in readreal:',
     -  icol1,' - ',icol2
      end if
      return
100   write (6,200) icol1,icol2,line
200   format(' ERROR: Invalid real was read between columns ',
     -  i3,' and ',i3,':',/,a132)
      return
1003  format(f3.0)
1004  format(f4.0)
1005  format(f5.0)
1006  format(f6.0)
1007  format(f7.0)
1008  format(f8.0)
1009  format(f9.0)
1010  format(f10.0)
1011  format(f11.0)
1012  format(f12.0)
1013  format(f13.0)
1014  format(f14.0)
1015  format(f15.0)
1016  format(f16.0)
1017  format(f17.0)
1018  format(f18.0)
1019  format(f19.0)
1020  format(f20.0)
      end
