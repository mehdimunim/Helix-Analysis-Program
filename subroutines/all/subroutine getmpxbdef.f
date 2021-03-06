      subroutine getmpxbdef(nslt,indexa,indexov,indexn,segid4,iresno,
     -   molsltlim,nsegslt,nanchorr,nanchorn,iout,maxrsd)
      dimension indexa(nslt),indexov(nslt),indexn(nslt),iresno(nslt),
     -  molsltlim(3,maxrsd)
      character*4 segid4(nslt)
      character*1 ansrun
c     print *,'GETMPXBDEF nslt,nsegslt=',nslt,nsegslt,' maxrsd=',maxrsd
      call zeroiti(indexa,0,nslt)
9134  print *,'Define the FIRST contact set'
      maxanchorlist=nslt
      ansrun=' '
      iall=0
      nanchorr=0
      do while (iall .eq. 0 .and. ansrun .ne. 'q')
        call quiz(ansrun,iansrun,' ',' ',1,
     -    'hydrophobic/salt bridge/contact anchor atoms',44,0,5,6,75)
        ifail=0
        call definelist(ansrun,nslt,nanchorr,indexa,indexn,nsegslt,
     -    segid4,iresno,molsltlim,'mpx contact',11,iall,maxanchorlist)
      end do
      print *,'Define the SECOND contact set'
      ansrun=' '
      iall=0
      nanchorn=0
      do while (iall .eq. 0 .and. ansrun .ne. 'q')
        call quiz(ansrun,iansrun,' ',' ',1,
     -    'hydrophobic/salt bridge/contact anchor atoms',44,0,5,6,75)
        ifail=0
        call definelist(ansrun,nslt,nanchorn,indexov,indexn,nsegslt,
     -    segid4,iresno,molsltlim,'mpx contact',11,iall,maxanchorlist)
      end do
c     Check for overlap
      nov=0
      do ia=1,nanchorr
        do ja=1,nanchorn
          if (indexa(ia) .eq. indexov(ja)) then
             if (nov .lt. 26) print *,'Atoms ', indexa(ia),' and ',
     -         indexov(ja),' are ','in both lists'
             nov=nov+1
          end if
        end do
      end do
      if (nov .gt. 0) then
        write (6,1002) nov
        call askstop(1)
        go to 9134
      end if
      write (6,1001) nanchorr,nanchorn
      write (iout,1001) nanchorr,nanchorn
      return
1001  format(' Number of reference atoms=',i6,/,
     -  ' Number of neighbour atoms=',i6)
1002  format(' Number of atoms in both lists=',i6,/,
     -  ' Mutually proximal contact list can not overlap')
      end
