      subroutine set_pdbinfo(iotyp,iwriteatsym,iwritecon,iobpdb,iocpdb,
     -  iaskcon)
      if (iotyp .eq. iobpdb .or. iotyp .eq. iocpdb) then
        if (iwriteatsym .lt. 0) call askyn(
     -    'Do you want chemical names written in the PDB file',50,
     -    1,-1,iwriteatsym,0,0)
        if (iaskcon .eq. 1) call askyn(
     -    'Do you want CONECT records written in the PDB file',50,1,
     -    -1,iwritecon,0,0)
      end if
      return
      end
