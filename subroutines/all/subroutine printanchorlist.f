      subroutine printanchorlist(label,llabel,ibondtyp,ianchor,nanchor,
     -  indexa,indexov,nanchorr,nanchorn,line,index,iresno,isegno,
     -  segid4,inamcol1,inamcol2,irescol1,irescol2,iout,mxrsd,mxrec)
      dimension ianchor(nanchor),indexa(mxrec),indexov(mxrec),
     -  index(mxrec),iresno(mxrec),isegno(mxrec),llabel(6)
      character*(*) label(6)
      character*4 segid4(mxrsd)
      character*132 line(mxrec)
c     print *,'PRINTANCHORLIST ibondtyp=',ibondtyp,' iout=',iout
      write (iout,*)
      if (ibondtyp .ne. 5) then
        write (iout,*) '=== List of ',
     -    label(ibondtyp)(1:llabel(ibondtyp)),' anchor atom indices:'
        call condenselist(ianchor,nanchor,0,iout)
        write (iout,*) '=== List of ',
     -    label(ibondtyp)(1:llabel(ibondtyp)),' anchor atoms:'
        do ia=1,nanchor
          ib=ianchor(ia)
          write (iout,1000) ia,line(index(ib))(irescol1:irescol2),
     -      iresno(ib),segid4(isegno(ib)),
     -      line(index(ib))(inamcol1:inamcol2),ib
        end do
      else
        write (iout,1001) 'reference',nanchorr
        call condenselist(indexa,nanchorr,0,iout)
        write (iout,1001) 'neighbour',nanchorn
        call condenselist(indexov,nanchorn,0,iout)
      end if
      write (iout,*)
      return
1000  format(i5,' Residue=',a,' (',i5,') Chain/seg=',a,' Atom=',a,
     -  ' (',i6,')')
1001  format(' Number of ',a,' atoms=',i6,' List:')
      end
