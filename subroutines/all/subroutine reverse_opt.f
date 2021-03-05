      subroutine reverse_opt(c,n,line,maxat)
      character*132 line(maxat)
      dimension c(3,maxat)
      dimension shift(3),rot(3,3)
      nlfound=0
      do il=1,100
        do icc=1,10
          if (line(il)(icc:icc+3) .eq. 'RTTR' .or.
     -        line(il)(icc:icc+3) .eq. 'TRRT') then
            nlfound=nlfound+1
            ic=icc+4
            if (line(il)(ic:ic) .eq. 't') then
              ic=ic+5
              read (line(il)(ic+1:ic+45),2000) shift
            else if (line(il)(ic:ic) .eq. 'r') then
              ic=ic+1
              read(line(il)(ic:ic),*) i
              ic=ic+5
              read (line(il)(ic+1:ic+33),2001) (rot(i,k),k=1,3)
            end if
          end if
        end do
c       print *,line(il)(1:77)
c       print *,'nlfound=',nlfound
      end do
      if (nlfound .lt. 4) go to 999
      write (6,1000) ((rot(i,k),k=1,3),i=1,3),shift
      call shiftmol(c,n,shift,c,1.0)
      call rotate_c(c,n,rot,c,'REVERSE',7)
      return
999   print *,'ERROR: no transformation information was found'
      return
1000  format(' Structure will be shifted by',5x,3f15.6,/,
     -  ' Structure will be rotated by',/,3(5x,3f11.7,/))
2000  format(3f15.6)
2001  format(3f11.7)
      end
