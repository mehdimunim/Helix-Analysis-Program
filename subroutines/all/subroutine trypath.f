      subroutine trypath(datapathtry,ldatapathtry,datapath,ldatapath)
      character*100 datapathtry,datapath
      character*200 filename
      filename(1:ldatapathtry)=datapathtry(1:ldatapathtry)
      filename(ldatapathtry+1:ldatapathtry+12)='/pdb_nam.dat'
      namlen=ldatapathtry+12
      call openfile(88,0,' ',1,'old',filename,namlen,notfound,2,1,0,1,
     -  0)
      if (notfound .eq. 0) then
        datapath=datapathtry
        ldatapath=ldatapathtry
        close (88)
        return
      else
        ldatapath=0
      end if
      return
      end
