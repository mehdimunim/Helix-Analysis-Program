      subroutine checknnlist(nfrst,n,ineig,nneig,nerr,maxng)
      dimension nneig(n),ineig(maxng,n)
c     Check for consistency of the neighbpor list
      nerr=0
      do ia=nfrst,n
        do ja=1,nneig(ia)
          ib=ineig(ja,ia)
          if (ib .lt. 1 .or. ib .gt. n) then
            print *,'PROGRAM ERROR: atom',ia,' has an invalid ',
     -        'bonded partner'
            nerr=nerr+1
          else if (nneig(ib) .eq. 0) then
            print *,'PROGRAM ERROR: atom',ib,' is bonded to',ia
            print *,'but it has no bonds listed'
            nerr=nerr+1
          else
            do jja=ja+1,nneig(ia)
              if (ineig(jja,ia) .eq. ib) then
                print *,'PROGRAM ERROR: ',
     -            'Duplicate bond partner of',ia,':',ib
                nerr=nerr+1
              end if
            end do
            do ii=1,nneig(ib)
              if (ineig(ii,ib) .eq. ia) go to 200
            end do
            print *,'PROGRAM ERROR: atom ',ia,' is not bonded to',ib
            nerr=nerr+1
200         continue
          end if
        end do
      end do
      if (nerr .gt. 0)
     -  print *,'NN list check for',n,' atoms done # of errors=',nerr
      return
      end
