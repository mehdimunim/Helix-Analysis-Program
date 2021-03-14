      subroutine set_hbondlims(hblimfac_def,hblimfac,angmin_def,
     -  angmin,iout)
      call getreal('H-bond length tolerance factor',30,hblimfac_def,
     -  hblimfac,1,23)
      call getreal('A-H...B angle minimum accepted',30,angmin_def,
     -  angmin,1,24)
      if (iout .gt. 0) write (iout,1000) hblimfac,angmin
      return
1000  format(' Hydrogen-bond is defined with hblimfac=',f5.2,
     -  ' and H-bond angle minimum=',f5.1)
      end
