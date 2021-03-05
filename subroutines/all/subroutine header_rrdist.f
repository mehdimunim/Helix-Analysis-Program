      subroutine header_rrdist(iw,nrefrange,nnegrange,resdistlim,
     -  irefseg1,ixrefres1,irefseg2,ixrefres2,irefres1,irefres2,
     -  inegseg1,ixnegres1,inegseg2,ixnegres2,inegres1,inegres2,
     -  incsolvrr,listrefres,listnegres,nrefres,nnegres,itemp,maxrsd)
      dimension listrefres(maxrsd),listnegres(maxrsd),itemp(maxrsd)
      write (iw,*)
      write (iw,2000) 'Reference',irefseg1,ixrefres1,irefseg2,
     -  ixrefres2,irefres1,irefres2
      if (nrefrange .gt. 1) then
        write (iw,2001) 'Reference'
        call zeroiti(itemp,0,maxrsd)
        do i=1,nrefres
          itemp(listrefres(i))=1
        end do
        call condensemask(itemp,1,listrefres(nrefres),iw,maxrsd)
      end if
      write (iw,2000) 'Neighbour',inegseg1,ixnegres1,inegseg2,
     -  ixnegres2,inegres1,inegres2
      if (nnegrange .gt. 1) then
        write (iw,2001) 'Neighbor'
        call zeroiti(itemp,0,maxrsd)
        do i=1,nnegres
          itemp(listnegres(i))=1
        end do
        call condensemask(itemp,1,listnegres(nnegres),iw,maxrsd)
      end if
      write (iw,2003)
      if (resdistlim .gt. 0.0) write (iw,2002) resdistlim
      if (incsolvrr .eq. 1) write (iw,2004) 'included'
      if (incsolvrr .eq. -1) write (iw,2004) 'excluded'
2000  format(1x,a,' residue range:',/,
     -  ' chain',i3,' res ',i5,' - ', 'chain',i3,' res ',i5,
     -  ' (residue serial #s',i5,' - ',i5,')')
2001  format(' NOTE: ',a,' residue range has gaps - actual ',
     -  'ranges:')
2002  format(' Residue pairs will be printed when the calculated ',
     -  'distance is <',f6.2,' A')
2003  format(" Residues not defined as 'reference' or 'neighbour' ",
     -  "residues are called 'other'")
2004  format(' Solvents are ',a,' in the additional neighbor list')
      return
      end
