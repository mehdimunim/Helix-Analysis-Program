      subroutine atnamtoatnum(atomnam,mmctype)
c     Determine atomic number from an atomname
      character*8 atomnam
      character*1 anam1(20)
      character* 132 line
      data anam1 /'H',3*' ','B','C','N','O','F',5*' ','P','S',
     -  3*' ','K'/
      line(1:5)=atomnam
      line(6:6)='X'
      ic=1
      call nextchar(line,ic,132)
      if (ic .le. 8) then
        do i=1,20
          if (atomnam(ic:ic) .eq. anam1(i)) then
             mmctype=i
             return
          end if
        end do
      end if
      mmctype=0
      return
      end
