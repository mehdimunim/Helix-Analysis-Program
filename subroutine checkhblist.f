      subroutine checkhblist(n,ineig,nhbneig,maxng)
      dimension nhbneig(n),ineig(maxng,n)
c     Check for consistency of the H-bond list
c     print *,'CHECKHBLIST n=',n
      do ia=1,n
        do ja=1,nhbneig(ia)
          ihb=ineig(maxng+1-ja,ia)
          if (nhbneig(ihb) .eq. 0) then
            print *,'PROGRAM ERROR: atom',ihb,' is on the HB list',
     -        ' of',ia
            print *,'but it has no HB partner listed'
          else
            do jja=1,ja-1
              if (ineig(maxng+1-jja,ia) .eq. ihb)
     -          print *,'Duplicate H-bond partner of',ia,':',ihb
            end do
            do ii=maxng+1-nhbneig(ihb),maxng
              if (ineig(ii,ihb) .eq. ia) go to 200
            end do
            print *,'PROGRAM ERROR: atom ',ia,' is not on the HB list',
     -        ' of',ihb
200         continue
          end if
        end do
      end do
      return
      end
