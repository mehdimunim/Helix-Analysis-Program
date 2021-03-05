      subroutine setmolres(ifres,ilres,isegno,molresflag,
     -  molsltlim,nrescol,irescol1,irescol2,resnames,nresslt,nmolslt,
     -  nsegslt,nmolsltnoion,minresflag,index,indexa,indexs,line,
     -  maxrsd,maxrec)
      dimension ifres(maxrec),ilres(maxrec),isegno(maxrec),
     -  molresflag(maxrsd),
     -  molsltlim(3,maxrsd),index(maxrec),indexa(maxrec),indexs(maxrec)
      character* 132 line(maxrec)
      character*8 resnames(maxrsd)
      character*8 ionresnam(100),molresnam(100)
c     Check for ions - make them separate molecules
      write (6,2117)
      call askyn('Do you have ions in this system',31,1,-1,ions,18,0)
      if (ions .gt. 0) call getnamelist(ionresnam,nrescol,nions,
     -  'Ion residue name',16,100)
      write (6,2118)
      call askyn('Do you have molecular residues in this system',
     - 45,1,-1,imolres,19,0)
      if (imolres .gt. 0) call getnamelist(molresnam,nrescol,
     -  nmolres,'Molecular residue name',22,100)
      if (ions+imolres .gt. 0) then
        nions_found=0
        nmolres_found=0
c       Gather residue names
        if (ions .gt. 0) call zeroiti(indexa,0,nions)
        if (imolres .gt. 0) call zeroiti(indexs,0,nmolres)
        minresflag=2
        do ir=1,nresslt
          resnames(ir)(1:nrescol)=
     -      line(index(ifres(ir)))(irescol1:irescol2)
          if (ions .gt. 0) then
            do irr=1,nions
              if (resnames(ir)(1:nrescol) .eq.
     -            ionresnam(irr)(1:nrescol)) then
                molresflag(ir)=2
                indexa(irr)=indexa(irr)+1
                nions_found=nions_found+1
              end if
            end do
          end if
          if (imolres .gt. 0) then
            do irr=1,nmolres
              if (molresflag(ir) .eq. 0 .and.
     -            resnames(ir)(1:nrescol) .eq.
     -            molresnam(irr)(1:nrescol))  then
                molresflag(ir)=1
                indexs(irr)=indexs(irr)+1
                nmolres_found=nmolres_found+1
              end if
            end do
          end if
          if (minresflag .gt. molresflag(ir)) minresflag=molresflag(ir)
        end do
        if (ions .gt. 0) then
          if (nions_found .eq. 0) write (6,2115) 'ions'
          if (nions_found .gt. 0) write (6,2114) 'ions',nions_found
        end if
        if (imolres .gt. 0) then
          if (nmolres_found .eq. 0) write (6,2115) 'molecular residues'
          if (nmolres_found .gt. 0)
     -      write (6,2114) 'molecular residues',nmolres_found
        end if
      else
        minresflag=0
      end if
      nrnoion=nresslt
      nmolion=0
      if (ions .gt. 0) then
c       Last block of solute residues may be ions
        write (6,2121) (ionresnam(irr)(1:nrescol),'ion',
     -    indexa(irr),irr=1,nions)
        ir=nresslt
        do while (ir .gt. 1 .and. molresflag(ir) .eq. 2)
          ir=ir-1
        end do
        if (ir .eq. 1 .and. molresflag(ir) .eq. 2) ir=0
        nrnoion=ir
        numresions=nresslt-nrnoion
        do ir=1,nrnoion
          if (molresflag(ir) .eq. 2) then
            write (6,2116) ir
            molresflag(ir)=1
          end if
        end do
        nsegnoion=isegno(ilres(nrnoion))
        do ir=nrnoion+1,nresslt
          molsltlim(1,nsegnoion+ir-nrnoion)=ifres(ir)
          molsltlim(2,nsegnoion+ir-nrnoion)=ilres(ir)
          molsltlim(3,nsegnoion+ir-nrnoion)=0
        end do
        nmolion=nresslt-nrnoion
        nmolslt=nmolslt+nmolion-(nsegslt-nsegnoion)
c        write (77,7272) (is,(molsltlim(k,is),k=1,2),is=1,nmolslt)
c7272    format(i4,' molsltlims=',2i6)
      end if
      if (imolres .gt. 0) then
        write (6,2121) (molresnam(irr)(1:nrescol),'molecular residue',
     -    indexs(irr),irr=1,nmolres)
c       Make requested residues molecules
        nmoladd=0
        nresslt=0
        do ir=nrnoion,1,-1
          if (molresflag(ir) .eq. 1) then
            nresslt=nresslt+1
            is=isegno(ilres(ir))
            if (ilres(ir) .eq. molsltlim(2,is) .and.
     -          ifres(ir) .eq. molsltlim(1,is)) then
c             Do nothing
c
            else if (ilres(ir) .eq. molsltlim(2,is)) then
c             Last residue of the segment
              nmoladd=nmoladd+1
              do iss=nmolslt+nmoladd,is+2,-1
                molsltlim(1,iss)=molsltlim(1,iss-1)
                molsltlim(2,iss)=molsltlim(2,iss-1)
              end do
              molsltlim(1,is+1)=ifres(ir)
              molsltlim(2,is+1)=ilres(ir)
              molsltlim(2,is)=ifres(ir)-1
            else if (ifres(ir) .eq. molsltlim(1,is)) then
c             First residue of the segment
              nmoladd=nmoladd+1
              do iss=nmolslt+nmoladd,is+2,-1
                molsltlim(1,iss)=molsltlim(1,iss-1)
                molsltlim(2,iss)=molsltlim(2,iss-1)
              end do
              molsltlim(1,is+1)=ilres(ir)+1
              molsltlim(2,is+1)=molsltlim(2,is)
              molsltlim(1,is)=ifres(ir)
              molsltlim(2,is)=ilres(ir)
            else
c             Residue in the middle of the segment
              nmoladd=nmoladd+2
              do iss=nmolslt+nmoladd,is+3,-1
                molsltlim(1,iss)=molsltlim(1,iss-2)
                molsltlim(2,iss)=molsltlim(2,iss-2)
              end do
              molsltlim(2,is)=ifres(ir)-1
              molsltlim(1,is+2)=ilres(ir)+1
              molsltlim(2,is+2)=molsltlim(2,is)
              molsltlim(1,is+1)=ifres(ir)
              molsltlim(2,is+1)=ilres(ir)
            end if
          end if
        end do
        nmolslt=nmolslt+nmoladd
c       write (77,7272) (is,(molsltlim(k,is),k=1,2),is=1,nmolslt)
      end if
      nmolsltnoion=nmolslt-nmolion
      if (ions+imolres .gt. 0) then
        do is=1,nmolslt
          molsltlim(3,is)=0
        end do
      end if
      return
2114  format(' Number of ',a,' found=',i4)
2115  format(' ERROR: none of the given ',a,' were found')
2116  format(' WARNING: residue ',i5,' is an ion residue but is not ',
     -  'at the end of the solute in a contiguous block ',/,
     -  ' - it will be treated as a molecular residue')
2117  format(' Simulaid treats ions as separate molecules for PBC ',
     -  'calculations',/,' and will not be considered part of the ',
     -  'solute ',/,' for the purpose of centering the system.',/,
     -  ' Ions have to be grouped together at the end of the solute ')
2118  format(' Simulaid treats segments as whole molecules for PBC ',
     -  'calculations.',/,' However, you may specify residue names ',/,
     -  ' that are to be considered separate molecules')
2121  format(' Number of ',a,1x,a,'s found=',i3)
      end
