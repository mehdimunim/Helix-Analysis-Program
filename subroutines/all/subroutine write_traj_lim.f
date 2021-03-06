      subroutine write_traj_lim(iout,label,llabel,i12,incr,iplot)
      character*(*) label
      character*200 trajnam12
      common /trajname/ trajnam12(2),ltrajnam12(2),ifirsttraj12(2),
     -  ilasttraj12(2),incrementtraj12(2)
      dimension llab12(2)
      character*8 lab12(2)
      data lab12 /'T','Second t'/,llab12 /1,8/
c     print *,'WRITE_TRAJ_NAME incr=',incr,' iplot=',iplot,' iout=',iout
      if (iplot .eq. 0) then
        if (llabel .gt. 0) write (iout,1000) label(1:llabel)
        write (iout,1001) lab12(i12)(1:llab12(i12)),
     -    trajnam12(i12)(1:ltrajnam12(i12)),
     -    ifirsttraj12(i12),ilasttraj12(i12),incrementtraj12(i12)
      else
        if (incrementtraj12(i12) .gt. 1) then
          write (iout,1002) lab12(i12)(1:llab12(i12)),
     -      trajnam12(i12)(1:ltrajnam12(i12)),
     -      ifirsttraj12(i12),ilasttraj12(i12),incrementtraj12(i12)
        else
          write (iout,1003) lab12(i12)(1:llab12(i12)),
     -      trajnam12(i12)(1:ltrajnam12(i12)),
     -      ifirsttraj12(i12),ilasttraj12(i12)
        end if
      end if
      incr=incrementtraj12(i12)
      return
1000  format(1x,a)
1001  format(1x,a,'rajectory analyzed: ',a,' Frames',i7,' - ',i6,
     -  ' Increment:',i5)
1002  format('(',a,'rajectory analyzed: ',a,' Frames',i7,' - ',i6,
     -  ' Increment:',i5,') show')
1003  format('(',a,'rajectory analyzed: ',a,' Frames',i7,' - ',i6,
     -  ') show')
      end
