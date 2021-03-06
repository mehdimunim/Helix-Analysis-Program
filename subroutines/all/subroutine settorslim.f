      subroutine settorslim(mask,molsltlim,nslt,nmolslt,iresno,segid4)
      dimension mask(nslt),molsltlim(3,nmolslt),iresno(nslt)
      character*4 segid4(nmolslt)
c     print *,'SETTORSLIM NSLT,NMOLSLT=',nslt,nmolslt
      call askyn('Do you want to set residue limits for torsions',46,
     -  1,-1,limres,0,0)
      if (limres .gt. 0) then
        call zeroiti(mask,0,nslt)
        if (nmolslt .gt. 1) write (6,2000) (i,segid4(i),i=1,nmolslt)
        idone=0
        do while (idone .eq. 0)
          imol=1
          if (nmolslt .gt. 1) then
            call getint('Molecule number (0 to finish)',29,imol,1,
     -        nmolslt,imol,000)
            if (imol .eq. 0) idone=1
          end if
          if (imol .gt. 0) then
            molf=molsltlim(1,imol)
            moll=molsltlim(2,imol)
100         call getrange(ifst,iresno(molf),ilst,iresno(moll),1,0,
     -       'residue to move (0 to finish)',29,iresno(moll),0)
            if (ifst*ilst .gt. 0) then
              iaf=molf
              do while (iresno(iaf) .lt. ifst .and. iaf .lt. moll)
                  iaf=iaf+1
              end do
              if (iresno(iaf) .ne. ifst) then
                write (6,*) 'ERROR: residue ID ',ifst,' is not found'
                go to 100
              end if
              ial=iaf
              do while (iresno(ial) .lt. ilst .and. ial .lt. moll)
                ial=ial+1
              end do
              if (iresno(ial) .ne. ilst) then
                write (6,*) 'ERROR: residue ID ',ilst,' is not found'
                go to 100
              end if
              do while (iresno(ial) .eq. ilst .and. ial .lt. moll)
                ial=ial+1
              end do
              if (iresno(ial) .ne. ilst) ial=ial-1
              write (6,2001) ifst,ilst,iaf,ial
              do ia=iaf,ial
                mask(ia)=1
              end do
            else
              idone=1
            end if
          end if
        end do
      else
        do ia=1,nslt
          mask(ia)=1
        end do
      end if
      return
2000  format(' Molecule/segment IDs:',(5(i3,1x,a4)))
2001  format(' Residue range [',i5,',',i5,'] includes the atom range [',
     -  i7,',',i7,']')
      end
