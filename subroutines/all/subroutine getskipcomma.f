      subroutine getskipcomma(inpt,line,len,llab,ifail)
      character*1000 line
      ifail=0
      lc=1
      do while (lc .le. 1)
        call blankout(line,1,len)
        read (inpt,1000,end=999) line
        call lastchar(line,lc,len)
      end do
c     print *,'GETSKIPCOMMA lc=',lc,' line:',line(1:lc)
      ndel=0
      do ic=1,lc
        if (line(ic:ic) .eq. ',' .or. line(ic:ic) .eq. ':') then
          ndel=ndel+1
        else
          line(ic-ndel:ic-ndel)=line(ic:ic)
        end if
      end do
      llab=lc-ndel
      return
1000  format(a)
999   ifail=1
      end
