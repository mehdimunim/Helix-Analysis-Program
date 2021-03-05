      subroutine asktrajform(inptrajtyp,ioutrajtyp,mmctrajtyp,
     -  resnamslv,inout,iasktrajtyp)
      character*8 resnamslv
      character*5 crdext
      common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
     -  iommod,iommc,iommc4,iogro,iomol2,iomae,iocif,ioxxx,ioins,
     -  ionxyz,iosxyz,iosxyzrq,iograsp,iofull,lext(19),crdext(19)
      character*8 resnamsv
      common /names/ resnamsv(19)
      character*1 ansfrm
      character*6 inoutlab
      data itrajtyp /0/
      mmctrajtyp=0
      if (iasktrajtyp .eq. 1) then
        inoutlab='output'
        if (inout .eq. -1) inoutlab=' input'
        call quiz(ansfrm,ians,' ',inoutlab,6,'trajectory file format',
     -    22,0,5,6,0)
        if (ansfrm .eq. 'c') then
          itrajtyp=1
        else if (ansfrm .eq. 'a') then
          itrajtyp=2
        else if (ansfrm .eq. 'l') then
          itrajtyp=3
        else if (ansfrm .eq. 'o') then
          itrajtyp=4
        else if (ansfrm .eq. 'x') then
          itrajtyp=5
        else if (ansfrm .eq. 'd') then
          itrajtyp=6
          print *,'Amber CDF trajectory input is not yet implemented'
          print *,'The Linux version of VMD can convert it to old Amber'
          stop
        end if
        if (inout .eq. +1) then
          ioutrajtyp=itrajtyp
        else
          inptrajtyp=itrajtyp
        end if
      end if
      if (itrajtyp .eq. 3 .or. iasktrajtyp .eq. 0) then
        call quiz(ansfrm,ians,' ',' ',0,'MMC trajectory file format',26,
     -    0,5,6,0)
        if (ansfrm .eq. 'b') then
          mmctrajtyp=1
        else if (ansfrm .eq. 'a') then
            mmctrajtyp=2
        else if (ansfrm .eq. 'n') then
          mmctrajtyp=3
        else if (ansfrm .eq. 'p') then
          mmctrajtyp=4
          resnamslv=resnamsv(iobpdb)
        else if (ansfrm .eq. 'c') then
          mmctrajtyp=5
          resnamslv=resnamsv(iocha)
        end if
      end if
      return
      end
