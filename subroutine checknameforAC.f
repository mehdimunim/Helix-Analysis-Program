      subroutine checknameforAC(nslt,line,index,inamcol1,inamcol2,
     -  maxrec)
      dimension index(maxrec)
      character* 132 line(maxrec)
      common /logging/ logfile,ipredict
c     Change atomnames using A for aromatic carbons to C
      naromc=0
      ichangeatoc=0
      do ia=1,nslt
        ifc=inamcol1
        call nextchar(line(index(ia)),ifc,inamcol2)
        if (line(index(ia))(ifc:ifc) .eq. 'A') then
          naromc=naromc+1
          if (naromc .eq. 1) then
            if (ipredict .eq. 0) then
              print *,'Atom names starting with A may stand for ',
     -          'aromatic carbons'
              call askyn('Do you want to change A*** atoms to C***',40,
     -          1,1,ichangeatoc,0,0)
            else
              write (6,1001)
              ichangeatoc=1
            end if
          end if
          if (ichangeatoc .eq. 1) line(index(ia))(ifc:ifc)='C'
        end if
      end do
      if (naromc .gt. 0) write (6,1000) naromc
      return
1000  format(' Number of atomnames changed from A*** to C***=',i4)
1001  format(' Atoms with name A*** (assumed to be aromatic carbons) ',
     -  'are renamed to C***',/,
     -  ' To avoid this, do not make the input predictable')
      end
