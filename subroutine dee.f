      subroutine dee(nframe,isg_i,ax1,ax2,cent1,cent2,halflen1,halflen2,
     -  distss,disthh,isg2ij)
      dimension ax1(3),ax2(3),cent1(3),cent2(3),s1(3),s2(3),h1(3),h2(3)
      parameter (MAXNHX=12,MAXNHX2=(MAXNHX*(MAXNHX-1)))
      character*2 ap_pa,in_ex
      common /signs/ itmem,normhx,isg2(MAXNHX2),memdir(MAXNHX),ap_pa(3),
     -  in_ex(2)
      call paramx(cent1,ax1,-halflen1,s1)
      call paramx(cent1,ax1,+halflen1,h1)
      if (nframe .eq. 1) then
c       Establish 2nd axis direction wrt the first
        call paramx(cent2,ax2,-halflen2,s2)
        call paramx(cent2,ax2,+halflen2,h2)
        dds1s2=dist2(s1,s2)
        dds1h2=dist2(s1,h2)
        ddh1h2=dist2(h1,h2)
        ddh1s2=dist2(h1,s2)
        if (dds1s2 .le. dds1h2 .and. ddh1h2 .le. ddh1s2) then
          isg2(isg_i)=1
        else if (dds1s2 .gt. dds1h2 .and. ddh1h2 .gt. ddh1s2) then
          isg2(isg_i)=-1
        else
          if (dds1s2+ddh1h2 .lt. dds1h2+ddh1s2) then
            isg2(isg_i)=1
          else
            isg2(isg_i)=-1
          end if
        end if
      end if
      call paramx(cent2,ax2,-isg2(isg_i)*halflen2,s2)
      call paramx(cent2,ax2,+isg2(isg_i)*halflen2,h2)
      distss=sqrt(dist2(s1,s2))
      disthh=sqrt(dist2(h1,h2))
      isg2ij=isg2(isg_i)
      return
      end
