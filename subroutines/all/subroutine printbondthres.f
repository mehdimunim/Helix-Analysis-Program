      subroutine printbondthres(ialist,nalist,ctfac,bondminfac,
     -  iatnm2,ramax,iout)
      dimension ialist(nalist)
      character*2 iatnm2(99)
      dimension ramax(99)
      write (iout,1003)
      write (iout,1002) (ialist(i),iatnm2(ialist(i)),i=1,nalist)
      do i=1,nalist
        do j=i,nalist
          call decidebondcut(ialist(i),ialist(j),rlim)
          write (iout,1000) iatnm2(ialist(i)),iatnm2(ialist(j)),
     -      sqrt(rlim),iatnm2(ialist(i)),ramax(ialist(i)),
     -      iatnm2(ialist(j)),ramax(ialist(j))
          if (ctfac .gt. 0.0)
     -      write (iout,1001) sqrt(rlim*ctfac),sqrt(rlim*bondminfac)
        end do
      end do
      return
1000  format(1x,a2,' - ',a2,' bondlength threshold=',f5.3,
     - 2(' ramax(',a2,')=',f5.3))
1001  format(9x,'clash threshold=',f5.3,
     -  ' minimum acceptable bondlength=',f5.3)
1002  format(' Atoms used: ',15(i3,1x,a2))
1003  format(' Bond threshold for heavy atoms: max(ramax(i),ramax(j))',
     -  /,' Bond threshold with hydrogen: 0.7*max(ramax(i),ramax(j))')
      end
