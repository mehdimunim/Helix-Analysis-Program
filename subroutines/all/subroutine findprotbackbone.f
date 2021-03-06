      subroutine findprotbackbone(line,index,iresno,ia1,inamcol1,
     -  ica,ic,in,icb,icg,icd,nslt,icaonly,maxrec)
      character* 132 line(maxrec)
      dimension index(maxrec),iresno(maxrec)
      character*4 atnami,atnamj
      ica=0
      ic=0
      in=0
      icb=0
      icg=0
      icd=0
      ia=ia1
      irespro=iresno(ia)
c     print *,'FINDPROTBACKBONE irespro,ia1=',irespro,ia1
      do while (iresno(ia) .eq. irespro .and. ia .le. nslt)
        atnami=line(index(ia))(inamcol1:inamcol1+3)
        call leftadjust4(atnami,atnamj)
        if (atnamj(1:2) .eq. 'CA') ica=ia
        if (atnamj(1:2) .eq. 'C ') ic=ia
        if (atnamj(1:2) .eq. 'N ') in=ia
        if (atnamj(1:2) .eq. 'CB') icb=ia
        if (atnamj(1:2) .eq. 'CG') icg=ia
        if (atnamj(1:2) .eq. 'CD') icd=ia
        ia=ia+1
      end do
c     print *,'ica,ic,in,icaonly=',ica,ic,in,icaonly
      if (icaonly .eq. -1) then
c       Full backbone is required
        if (ica*ic*in .eq. 0) then
          print *,'ERROR: residue ',irespro,' has no CA or C or N'
          stop
        end if
      else if (icaonly .eq. 1) then
c       Only CA is required
        if (ica .eq. 0) then
          print *,'ERROR: residue ',irespro,' has no CA'
          stop
        end if
      else
        if (ica .gt. 0 .and. ic*in .eq. 0) then
c         N and/or C is missing, CA is present - ask for CA only option
          print *,'Residue ',irespro,' has no C or N but CA is present'
          call askyn('Do you want to look for just CAs',32,1,+1,icaonly,
     -      0,0)
          if (icaonly .eq. 0) stop
        end if
      end if
      if (ic .eq. 0) ic=ica
      if (in .eq. 0) in=ica
      ia1=ia
      return
      end
