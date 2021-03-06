      subroutine writeuniquelist(index,ixresno,n,resnames,nrcol,iout,
     -  ires1,ires2,iresinc,itemp1,itemp2,label,llabel,mxrsd)
      dimension index(mxrsd),ixresno(mxrsd),itemp1(mxrsd),itemp2(mxrsd)
      character*(*) label
      character*8 resnames(mxrsd)
c     print *,'WRITEUNIQUELIST ires1,ires2=',ires1,ires2,
c    -  ' lab=',label(1:llabel)
      nmem=0
      call zeroiti(itemp1,0,n)
      call zeroiti(itemp2,0,n)
      do i=ires1,ires2
        if (index(i) .ne. 0) then
          nmem=nmem+1
          itemp1(nmem)=i
        else
          itemp2(i)=i
        end if
      end do
c     if (nmem .eq. 0) print *,'WRITEUNIQUELIST nmem=',nmem,
c    -  ' ires1,ires2=',ires1,ires2
      if (nmem .eq. 0) then
        write (iout,1003) 'contact',label(1:llabel)
      else
        write (iout,1000) 'contact',label(1:llabel),
     -   (' (',resnames(itemp1(i))(1:nrcol),ixresno(itemp1(i)),i=1,nmem)
        write (iout,1001) 'contact',label(1:llabel)
        call condensemask(index,ires1,ires2,iout,mxrsd)
        if (iresinc .gt. 0) write (iout,1002) iresinc
      end if
      nnonmem=0
      do i=1,n
        if (itemp2(i) .ne. 0) then
          nnonmem=nnonmem+1
          itemp1(nnonmem)=i
        end if
      end do
      if (nnonmem .gt. 0) then
        write (iout,1000) 'non-contact',label(1:llabel),
     -    (' (',resnames(itemp1(i))(1:nrcol),ixresno(itemp1(i)),
     -    i=1,nnonmem)
        write (iout,1001) 'non-contact',label(1:llabel)
        call condensemask(itemp2,ires1,ires2,iout,mxrsd)
        if (iresinc .gt. 0) write (iout,1002) iresinc
      end if
      return
1000  format(' Nonredundant list of ',a,1x,a,' residues:',/,
     -  5(a,a,i5,')'))
1001  format(' Residue SEQUENCE numbers of the ',a,1x,a,' residues:')
1002  format(' Residue numbers are decremented by',i5,' to give the ',
     -  'sequence numbers')
1003  format(' List of ',a,1x,a,' residue numbers is empty')
      end
