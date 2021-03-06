      subroutine setcolor(ic)
      common /depthcuedat/ near,ifar,ramp0,idepth,idepthon,idrawh,
     -  linew,isidec,nbackb,idrawslv
c     print *,'SETCOLOR idepthon,ic,idepth=',idepthon,ic,idepth
      if (idepthon*idepth .eq. 1) then
        if (ic .eq. 0) then
        else
          i1=100*ic+1
          i2=100*(ic+1)
        end if
      else
      end if
      return
      end
