      subroutine setrepats(nsegslt,molsltlim,nslt,nneig,ineig,
     -  molresflag,minresflag,it1,it2,maxneig,maxrsd)
      dimension molsltlim(3,maxrsd),nneig(nslt),ineig(maxneig,nslt),
     -  it1(nslt),it2(nslt),molresflag(maxrsd)
      dimension moltyp(3)
      character*1 ctyp
      character*41 question
c     print *,'SETREPATS  maxneig,maxrsd=', maxneig,maxrsd
c     Establish solute molecule centers (if required)
      call zeroiti(moltyp,0,3)
      do is=1,nsegslt
        moltyp(molresflag(is)+1)=moltyp(molresflag(is)+1)+1
      end do
      if (moltyp(3) .gt. 0 .and. minresflag .eq. 0) then
        call askyn(
     - 'Do you want to include the molecular residues in the centering',
     -    62,1,1,molrescent,80,0)
        if (molrescent .eq. 1) minresflag=minresflag+1
      end if
c     molsltlim(1,im), molslt(2,im): first and last atom of solute molecule im
c     molsltlim(3,im): atom number of topological center;
c                    zero if geom center was selected; -1 if COM is selected
      if (nsegslt .le. 1) then
        nsegslt=1
        molsltlim(1,1)=1
        molsltlim(2,1)=nslt
        molsltlim(3,1)=0
      end if
c     Optionally establish representative atoms
      call quiz(ctyp,ictyp,'g',' ',0,'molecular center',16,0,5,6,25)
      if (ctyp .eq. 'g') then
        do is=1,nsegslt
          molsltlim(3,is)=0
        end do
      else if (ctyp .eq. 'c') then
        do is=1,nsegslt
          molsltlim(3,is)=-1
        end do
      else
        iout=6
        do is=1,nsegslt
          if (ctyp .eq. 't' .and. is .le. 5) iout=6
          call findtcent(ineig,nneig,it1,it2,molsltlim(1,is),
     -      molsltlim(2,is),molsltlim(3,is),iout,maxneig,nslt)
          if (ctyp .eq. 't' .and. is .eq. 6) then
            print *,'Topology center list truncated'
            iout=0
          end if
        end do
        if (ctyp .eq. 'i') then
          do is=1,nsegslt
            write (question,1000) is
100         call getint(question,41,molsltlim(3,is),1,molsltlim(2,is),
     -        molsltlim(3,is),0)
            if (molsltlim(3,is) .lt. molsltlim(1,is) .or.
     -          molsltlim(3,is) .gt. molsltlim(2,is)) then
              print *,'Invalid number'
              go to 100
            end if
          end do
        end if
      end if
      return
1000  format('PBC center atom for solute molecule #',i4)
      end
