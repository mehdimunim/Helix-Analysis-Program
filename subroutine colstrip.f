      subroutine colstrip(cvav,nres,ixdel,scalefac,iydown,nrep,lab,ncol,
     -  ips)
      dimension cvav(nres)
      character*2 lab
      riy=iydown
      if (nrep .le. 1) then
        write (ips,3004) 12
        write (ips,*) 'np'
        write (ips,3001) ixdel,-iydown
      end if
      kcprev=0
      ixcurr=kcprev
      ix=1
      rrprev=0.0
      ixxprev=ixdel
      do while (ix .lt. nres)
        do while (ixcurr .eq. kcprev .and. ix .lt. nres)
          ix=ix+1
          ixcurr=ncol*cvav(ix)+1
          if (ix .eq. nres) ixcurr=0
        end do
        if (kcprev .gt. ncol) kcprevncol=ncol
        rr=ix-1
        if (nrep .le. 1) then
          call rgbcolor(ips,kcprev)
          write (ips,*) 'np'
          ixx=ixdel+scalefac*(ix-1)
          write (ips,3001) ixxprev,-iydown
          write (ips,3002) ixx,-iydown
          write (ips,*) 'sk'
          ixxprev=ixx
        end if
        kcprev=ixcurr
      end do
      rr=nres
          rx=nres+5
      if (nrep .le. 1) then
        call rgbcolor(ips,kcprev)
        write (ips,*) 'np'
        ixnres=ixdel+scalefac*nres
        write (ips,3001) ixnres,-iydown
        call rgbcolor(ips,9)
        write (ips,3001) ixnres+5,-iydown-5
        call psshow(ips,lab,2)
        write (ips,*) 'sk'
      end if
      return
3001  format(i4,i5,' m')
3002  format(i4,i5,' l')
3004  format(i5,' setlinewidth')
      end
