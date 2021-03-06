      subroutine molarvol(nslt,n,naslv,line,index,ir1,ir2,irn1,irn2,
     -  naa,nna,nnw,nnf,v,maxrec)
c     Calculate molar volume using known values of amino and nucleic acids
      character* 132 line(maxrec)
      dimension index(maxrec)
      character*8 resname,resnamela
      character*1 aanames1
      character*2 mmodtoamb
      character*3 aanames3
      common /atnamcon/ mmodtoamb(100),aanames1(58),aanames3(58),
     -  naanames,nnanames,nnammnames,nnames,ixwatnam
      common /pmvol/ pmvaana(58)
      lresn=ir2-ir1+1
c     print *,'MOLARV ir1,ir2,irn1,irn2=',ir1,ir2,irn1,irn2
c     print *,'MOLARV nslt=',nslt
      vna=0.0
      vpr=0.0
      vw=0.0
      naa=0
      nna=0
      nnw=0
      nnf=0
      iresnumo=0
      do ia=1,nslt
        call readint(line(index(ia)),irn1,irn2,iresnum,2,1,irerr)
        if (iresnum .ne. iresnumo) then
c         New residue starts
          iresnumo=iresnum
          resname(1:lresn)=line(index(ia))(ir1:ir2)
          call leftadjustn(resname,resnamela,lresn)
          call findname(resnamela,aanames3,1,naanames,ix,3)
          if (ix .gt. 0 .and.
     -        ix .ne. ixwatnam-1 .and. ix .ne. ixwatnam) then
            vpr=vpr+pmvaana(ix)
            naa=naa+1
          else if (resname(1:3) .eq. 'TIP' .or.
     -       resname(1:3) .eq. 'WTR' .or. resname(1:3) .eq. 'HOH') then
            vw=vw+pmvaana(32)
            nnw=nnw+1
          else
            call findname(resnamela,aanames3,
     -        naanames+nnammnames+1,nnames,ix,3)
            if (ix .gt. 0) then
              vna=vna+pmvaana(ix)
              nna=nna+1
            else
              nnf=nnf+1
            end if
          end if
c         write (77,*) 'ia=',ia,' resname=',resname(1:3),
c    -      ' naa,nnw,nna,nnf=',naa,nnw,nna,nnf
        end if
      end do
      if (n .gt. nslt) then
        resname(1:lresn)=line(index(nslt+1))(ir1:ir2)
        call leftadjustn(resname,resnamela,lresn)
        if (resnamela(1:3) .eq. 'TIP' .or.
     -      resname(1:3) .eq. 'WTR' .or. resname(1:3) .eq. 'HOH') then
          nnw=(n-nslt)/naslv
          vw=nnw*18.07*1.e+24/6.022045e+23
        end if
      end if
      if (nslt .gt. 0) then
        vpr=vpr*1.e+24/6.022045e+23
        vna=vna*1.e+24/6.022045e+23
        v=vpr+vna
        write (6,1000) naa,nna,nnf,v
        if (naa .gt. 0) write (6,1001) 'protein',vpr
        if (nna .gt. 0) then
          write (6,1001) 'nucleic acid',vna
          write (6,1003)
        end if
        if (nnf .gt. 0) write (6,1002)
      end if
      if (n .gt. nslt) then
        if (nnw .gt. 0) then
          print *,'The solvent is assumed to be water'
          write (6,1001) 'water',vw
        end if
      end if
      return
1000  format(' The solute contains ',i4,' amino acid residues',
     -  i4,' nucleic acid residues',/,
     -  ' and',i6,' unclassified residues',/,
     -  ' The volume of the solute is estimated to be ',f10.2,' A^3')
1001  format(' Volume of the ',a,' (part) is estimated to be',f12.2,
     -  ' A^3')
1002  format(' NOTE: unclassified residues did not contribute to the ',
     -  'total solute volume')
1003  format(' NOTE: Nucleic acid volumes are very sensitive to ',
     -  'cc, salt, etc.')
      end
