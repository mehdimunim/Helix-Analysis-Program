      subroutine checkdupnames(namelist,lnamelist,n,label,llabel,idup)
      character*(*) namelist,label
      dimension namelist(n)
c     Checks list of names for duplicates
      idup=0
      do i=2,n
        call findname(namelist(i),namelist,1,i-1,ix,lnamelist)
        if (ix .gt. 0) then
          idup=1
          if (llabel .gt. 0) write (6,1000) label(1:llabel),
     -        i,ix,namelist(i)(1:lnamelist)
        end if
      end do
      return
1000  format(' Duplicate found in list ',a,
     -  ': list(',i4,')=list(',i4,')=',a)
      end
