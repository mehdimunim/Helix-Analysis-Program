      subroutine askrot(nstep)
      common /rotmat/ matrot0(4,4),matrot(4,4),nomat0
      character*1 ans,ans0,anslc,xyzorig(3)
      common /askrotdat/ angle,ans0
      common /depthcuedat/ near,ifar,ramp0,idepth,idepthon,idrawh,
     -  linew,isidec,nbackb,idrawslv
      character*1 xyz
      common /axislab/ xyz(3)
      character*3 onoff(2)
      data onoff /'OFF','ON '/,xyzorig /'x','y','z'/
      if (nstep .eq. 0) then
        ans0='x'
        angle=10.0
      end if
      if (mod(nstep,20) .eq. 0)
     -  write (6,2104) onoff(2-idepth),onoff(2-idrawh),onoff(2-isidec),
     -  onoff(2-idrawslv),linew
      nstep=nstep+1
      ans=' '
100   read (5,1010) ans
      if (ans .eq. 's' .or. ans .eq. 'S') then
        nstep=-1
        return
      else if (ans .eq. 'q' .or. ans .eq. 'Q') then
        call getreal('New angle increment',19,angle,angle,0,36)
        go to 100
      else if (ans .eq. 'd' .or. ans .eq. 'D') then
        idepth=1-idepth
        write (6,2001) 'Depth cueing',onoff(idepth+1)
        go to 100
      else if (ans .eq. 'h' .or. ans .eq. 'H') then
        idrawh=1-idrawh
        write (6,2001) 'Hydrogens',onoff(idrawh+1)
        go to 100
      else if (ans .eq. 'c' .or. ans .eq. 'C') then
        if (nbackb .gt. 0) then
          isidec=1-isidec
          write (6,2001) 'Side chains',onoff(isidec+1)
        else
          print *,'There are no backbone atoms - input ignored'
        end if
        go to 100
      else if (ans .eq. 'v' .or. ans .eq. 'V') then
        idrawslv=1-idrawslv
        write (6,2001) 'Solvents',onoff(idrawslv+1)
        go to 100
      else if (idigit(ans,1) .eq. 1) then
        read (ans,*) linew
        write (6,2002) linew
      end if
      if (ans .eq. 'x' .or. ans .eq. 'y' .or. ans .eq. 'z'
     -    .or. ans .eq. 'X' .or. ans .eq. 'Y' .or. ans .eq. 'Z') then
        anslc=ans
        call uplow(ans,anslc,1,noabc)
        do k=1,3
          if (anslc .eq. xyz(k)) ans=xyzorig(k)
        end do
        ans0=ans
      else
        ans=ans0
      end if
      iaxis=10.0*angle
c     print *,'iaxis=',iaxis,' ans=',ans,' angle=',angle
      return
1010  format(a1)
2001  format(1x,a,' has been turned ',a)
2002  format(' Solute linewidth has been changed to',i2)
2104  format(' Type x/y/z to rotate by 10 degrees around axis x/y/z',/,
     -  ' Type q to change the angle increment',/,
     -  ' Type d to toggle depth-cueing ',a,/,
     -  ' Type h to toggle showing hydrogens  ',a,/,
     -  ' Type c to toggle showing sidechains ',a,/,
     -  ' Type v to toggle showing solvents   ',a,/,
     -  ' Type a number (1 digit) to change the linewidth from',i2,/,
     -  ' Hit <Enter> to continue rotating around the previous axis',/,
     -  ' Type s to stop the rotation and continue the run')
      end
