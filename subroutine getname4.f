      subroutine getname4(ic,line,name,len,lead_trail)
      character*(*) line
      character*4 name
      name='    '
      call nextchar(line,ic,len)
      icf=ic
      call nextblank(line,ic,len)
      if (ic-icf .le. 4 .or. lead_trail .eq. 1) then
        name(1:ic-icf)=line(icf:ic-1)
      else
c       Use trailing part of the name
        name=line(ic-4:ic-1)
      end if
      return
      end
