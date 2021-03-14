      subroutine writefree(line,ioutyp,chemnam,cx,cy,cz,ires,q)
      character*4 chemnam
      character* 132 line
      character*5 crdext
      common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
     -  iommod,iommc,iommc4,iogro,iomol2,iomae,iocif,ioxxx,ioins,
     -  ionxyz,iosxyz,iosxyzrq,iograsp,iofull,lext(19),crdext(19)
c     Limited free-format outputs for InsightII
      if (ioutyp .eq. ionxyz) then
        write (line,1100) ianum(chemnam,0,4),cx,cy,cz
      else if (ioutyp .eq. iosxyz) then
        write (line,1101) chemnam,cx,cy,cz
      else if (ioutyp .eq. iosxyzrq) then
        write (line,1102) chemnam,cx,cy,cz,ires,q
      end if
      return
1100  format(i5,3f10.5)
1101  format(a5,3f10.5)
1102  format(a5,3f10.5,i5,f10.5)
      end
