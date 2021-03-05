      subroutine PDBtommc(resn,atomn,ptype,char,indx,gcent,igr,
     -  ideoxymmc,ifile,nunknown,nunknownr)
      character*1 gcent,ans1
      character*4 ptype
      character*8 resn,resnam,atomn,atomnam
      character*4 ires
      character*8 convdat,grpinfo
      character*200 sfilename
      common /savedat/ mxresdat,maxcondat,ifst(1000),ilst(1000),
     -  nres,iresgen,lsfilename,ires(1000),convdat(7,10000),sfilename
c     print *,'PDBtommc resn=',resn,' atomn=',atomn
      ptype='****'
      char=0.0
      gcent=' '
      igr=0
      call leftadjustn(atomn,atomnam,8)
      call leftadjustn(resn,resnam,8)
c     If residue name/atom name not found, check alternatives
      if (ifile .eq. 2) then
c       Original Amber
        if (resnam .eq. 'ADE     ' .or. resnam .eq. 'GUA     ' .or.
     -    resnam .eq. 'CYT     ') then
          if (ideoxymmc .eq. -1) then
c           Ask if oxy or deoxy NA
            call getname(ans1,len,
     -        'Nucleic acid type (D for Deoxy, R Oxy)',38,1,'',0,0,0,0)
            ideoxymmc=1
            if (ans1 .eq. 'r' .or. ans1 .eq. 'R') ideoxymmc=0
          end if
          resnam(2:4)=resnam(1:3)
          if (ideoxymmc .eq. 0) resnam(1:1)='R'
          if (ideoxymmc .eq. 1) resnam(1:1)='D'
        else if (resnam .eq. 'THY  ') then
          ideoxymmc=1
          resnam='DTHY    '
        else if (resnam .eq. 'URA     ') then
          resnam='RURA    '
          ideoxymmc=0
        end if
      else if (ifile .eq. 3) then
c       Amber 94
      else if (ifile .eq. 4) then
c       Charmm
      end if
c???
c     if (idigit(atomnam(1:1),1) .eq. 1) then
c       atomnam(1:4)=atomnam(2:5)
c       atomnam(5:5)=' '
c     end if
      do i=1,nres
c       write (77,*) 'resnam,ires(i)=',resnam,ires(i),' if,l=',ifst(i),ilst(i)
        if (resnam(1:4) .eq. ires(i)) then
          do j=ifst(i),ilst(i)
c            write (77,1077)
c     -        resnam,atomnam,convdat(1,j)(1:4),convdat(2,j)(1:4)
c1077        format(' rn,an=',a4,',',a4,'*',' cd1,2=',a4,',',a4,'*')
            if (atomnam(1:4) .eq. convdat(2,j)(1:4)) then
              ptype=convdat(3,j)(1:4)
              read (convdat(4,j),1001) char
              grpinfo=convdat(5,j)
              go to 100
            end if
          end do
c         If not found, try ACE, NME or DPOM or RPOM
c         since Amber uses these groups as different residues
          do ii=1,nres
            if (ires(ii) .eq. 'ACE ' .or. ires(ii) .eq. 'NME ' .or.
     -          (ires(ii)(2:4) .eq. 'POM' .and.
     -           ires(ii)(1:1) .eq. resnam(1:1))) then
              do j=ifst(ii),ilst(ii)
                if (atomnam(1:4) .eq. convdat(2,j)(1:4)) then
                  ptype=convdat(3,j)(1:4)
                  read (convdat(4,j),1001) char
                  grpinfo=convdat(5,j)
                  go to 100
                end if
              end do
            end if
          end do
          if (iresgen .gt. 0) then
c           Check if there are residue-name independent atomnames
            do j=ifst(iresgen),ilst(iresgen)
              if (atomnam(1:4) .eq. convdat(2,j)(1:4)) then
                ptype=convdat(3,j)(1:4)
                read (convdat(4,j),1001) char
                grpinfo=convdat(5,j)
                go to 100
              end if
            end do
          end if
          write (6,2000) atomnam,resnam,indx
          nunknown=nunknown+1
          return
100       if (grpinfo(1:5) .eq. 'GROUP') then
            read(grpinfo(6:6),1003) igcnt
            if (igcnt .eq. 1) gcent='G'
            read(grpinfo(7:8),1004) igr
          end if
          return
        end if
      end do
      write (6,2001) resnam
      nunknownr=nunknownr+1
      return
1001  format(f8.0)
1003  format(i1)
1004  format(i2)
2000  format(' ERROR: atomname=',a5,' resname=',a5,' (',i5,')',
     -  ' not found in conversion table')
2001  format(' ERROR: resname=',a4,' not found in conversion ',
     -  'table')
      end
