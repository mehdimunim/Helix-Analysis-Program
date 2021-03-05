      subroutine hydropathylist(n,nslt,ixres,resnames,cv,ihydtyp,
     -  maxrsd,maxrec)
      character*8 resnames
      dimension ixres(maxrec),resnames(maxrsd),cv(n)
      character*3 hphres
      character*1 ans
      dimension hph3(4,23),hphres(23)
c     EIS KD WHI
      data hphres /
     -  'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS',
     -  'ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR',
     -  'VAL','ASX','GLX','UNK'/
      data hph3 /
     -  0.62,1.80,-0.50,0.0, -2.53,-4.50,-1.81,0.0,
     -  -0.78,-3.50,-0.85,0.0, -0.90,-3.50,-3.64,0.0,
     -  0.29,2.50,0.02,0.0, -0.85,-3.50,-0.77,0.0,
     -  -0.74,-3.50,-3.63,0.0, 0.48,-0.40,-1.15,0.0,
     -  -0.40,-3.20,-2.33,0.0, 1.38,4.50,1.12,0.0,
     -  1.06,3.80,1.25,0.0, -1.50,-3.90,-2.80,0.0,
     -  0.64,1.90,0.67,0.0, 1.19,2.80,1.71,0.0,
     -  0.12,-1.60,-0.14,0.0, -0.18,-0.80,-0.46,0.0,
     -  -0.05,-0.70,-0.25,0.0, 0.81,-0.90,2.09,0.0,
     -  0.26,-1.30,0.71,0.0, 1.08,4.20,0.46,0.0,
     -  -0.84,-3.50,-2.25,0.0, -0.80,-3.50,-2.20,0.0,
     -  0.00,-0.49,-0.52,0.0/
      if (ihydtyp .eq. 0) then
        call quiz(ans,ihydtyp,' ',' ',0,
     -  'hydropathy/hydrophobicity scale source',38,0,5,6,0)
        if (ihydtyp .eq. 4) then
          do i=1,23
            call getreal('Residue '//hphres(i)//' hydrophobicity',26,
     -        999999.0,hph3(4,i),0,0)
          end do
        end if
        write (6,1002) (hphres(i),(hph3(j,i),j=1,4),i=1,23)
      end if
      do ia=1,nslt
        ir=1
        do while (ir .lt. 23 .and.
     -     resnames(ixres(ia))(1:3) .ne. hphres(ir))
          ir=ir+1
        end do
        cv(ia)=hph3(ihydtyp,ir)
      end do
c     Set solvent atoms to unknown
      do ia=nslt+1,n
        cv(ia)=hph3(ihydtyp,23)
      end do
      return
1002  format(' The residue values stored:',/,
     -   '    Kyte-Dooli Eisenberg     White    Input',/,
     -  (1x,a3,4f10.2))
      end
