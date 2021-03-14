      subroutine getanchormod(ianchor2,iselfanc,nosameseg,iallsel,iout)
      character*8 inout(2)
      character*9 onetwo(2)
      data inout /'excluded','included'/,onetwo/'one end  ','both ends'/
      if (iallsel .eq. 0) then
        call askyn(
     -    'Do you want anchor atoms at both ends',37,1,-1,ianchor2,71,0)
        iselfanc=1
        if (ianchor2 .eq. 0)
     -    call askyn('Do you want exclude anchor-anchor bonds',39,0,-1,
     -      iselfanc,000,0)
      else
        ianchor2=0
        iselfanc=1
      end if
      call askyn('Do you want exclude intra segment/chain bonds',45,1,
     -  -1,nosameseg,000,0)
      if (iout .gt. 0) then
        write (iout,1000) 'Anchor-anchor',inout(iselfanc+1)
        write (iout,1000) 'Intra segment/chain',inout(2-nosameseg)
        write (iout,1001) onetwo(ianchor2+1)
      end if
      return
1000  format(1x,a,' bonds will be ',a)
1001  format(' Anchor atom is required at ',a,' of a bond')
      end
