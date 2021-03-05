      subroutine connlinix(p,i1,i2,ixdup,iconntyp,ioutpdb,max)
      dimension p(3,max)
      dimension ixdup(max)
c     Just draw a line connecting p(i1) & p(i2)
c     iconntyp = 1,2,3: SGI, PDB, both
      if (iconntyp .ne. 2) then
      end if
      if (iconntyp .ne. 1) then
        if (ixdup(i1) .lt. ixdup(i2)) then
          write (ioutpdb,1000) ixdup(i1),ixdup(i2)
        else
          write (ioutpdb,1000) ixdup(i2),ixdup(i1)
        end if
      end if
      return
1000  format('CONECT',2i5)
      end
