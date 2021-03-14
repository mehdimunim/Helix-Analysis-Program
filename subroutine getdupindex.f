      subroutine getdupindex(nsegm,segid4,index)
      character*4 segid4(nsegm)
      dimension index(nsegm)
c     print *,'GETDUPINDEX nsegm=',nsegm
      call indexit(index,1,nsegm,0)
      do is=2,nsegm
        do iss=1,is-1
          if (segid4(is) .eq. segid4(iss)) index(is)=iss
        end do
      end do
      return
      end
