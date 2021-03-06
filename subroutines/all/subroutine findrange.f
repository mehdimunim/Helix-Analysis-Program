      subroutine findrange(index,n1,n,ifind,ifirst,ilast,label,llabel,
     -  isilent,ifail)
c*****Find the atom range in index containig ifind
      dimension index(n)
      character*(*) label
c     print *,'FINDRANGE n1,n,ifind=',n1,n,ifind
      ifail=0
      ifirst=n1
      do while (index(ifirst) .ne. ifind)
        ifirst=ifirst+1
        if (ifirst .gt. n) then
          if (isilent .eq. 0) print *,'ERROR: subroutine findrange ',
     -      'failed to find ',label(1:llabel),' ',ifind
          ifail=1
          return
        end if
      end do
      ilast=ifirst
      do while (index(ilast) .eq. ifind .and. ilast .lt. n)
        ilast=ilast+1
      end do
      if (ilast .lt. n) ilast=ilast-1
c     print *,'FINDRANGE ifirst,ilast=',ifirst,ilast
      return
      end
