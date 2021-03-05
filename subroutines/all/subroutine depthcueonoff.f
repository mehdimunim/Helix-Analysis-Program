      subroutine depthcueonoff(ion)
      common /depthcuedat/ near,ifar,ramp0,idepth,idepthon,idrawh,
     -  linew,isidec,nbackb,idrawslv
      if (ion .eq. 1) then
        idepthon=1
      else
        idepthon=0
      end if
      return
      end
