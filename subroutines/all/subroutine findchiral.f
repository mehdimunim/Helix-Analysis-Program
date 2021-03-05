      subroutine findchiral(n,ian,nneig,nhneig,ineig,ichiral,maxneig)
      dimension ian(n),nneig(n),nhneig(n),ineig(maxneig,n),
     -  ichiral(n)
      dimension nnpairs(2,6)
      call zeroiti(ichiral,0,n)
      if (n .gt. 1000) then
        print *,'Chirality search is done for max 1000 atoms'
        return
      end if
      do ia=1,n
        if (nneig(ia) .eq. 4 .and. nhneig(ia) .le. 1) then
          call listsame(ineig(1,ia),ian,4,n,nnpairs,nsamepair)
          if (nsamepair .eq. 0) then
            ichiral(ia)=1
          else
            ipair=nsamepair
            ndiff=0
            do while (ipair .gt. 0)
              call comparetree(nneig,ineig,ian,ia,nnpairs(1,ipair),
     -          nnpairs(2,ipair),idiff,n,maxneig)
              if (idiff .eq. 0) then
                ipair=0
              else
                ipair=ipair-1
                ndiff=ndiff+1
              end if
            end do
            if (ndiff .eq. nsamepair) ichiral(ia)=1
          end if
        end if
      end do
      return
      end
