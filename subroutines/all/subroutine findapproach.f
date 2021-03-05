      subroutine findapproach(c,irra1,irra2,inra1,inra2,iatnum,ignoreh,
     -  irarepm,itarepm,rmin2,maxrec)
      dimension c(3,maxrec),iatnum(maxrec)
      rmin2=100000.0
      if (ignoreh .eq. 0) then
c       Use all atoms
        do ira=irra1,irra2
          do ita=inra1,inra2
            r2=dist2(c(1,ira),c(1,ita))
            if (r2 .lt. rmin2) then
              rmin2=r2
              irarepm=ira
              itarepm=ita
            end if
          end do
        end do
      else
c       Use the heavy atoms only
        do ira=irra1,irra2
          if (iatnum(ira) .gt. 1) then
            do ita=inra1,inra2
              if (iatnum(ita) .gt. 1) then
                r2=dist2(c(1,ira),c(1,ita))
                if (r2 .lt. rmin2) then
                  rmin2=r2
                  irarepm=ira
                  itarepm=ita
                end if
              end if
            end do
          end if
        end do
      end if
      return
      end
