      subroutine setpbccell(quest,lquest,edge,edge_gen,cell,ncell,
     -  cellalt,ixyzhex,npbc,ioppbc,iusepbc,vol,nw,rinscr,
     -  rcirc,nonone)
      dimension edge(3),edge_gen(3,3),cell(3,27),cellalt(3,27),
     -  ixyzhex(3)
      character*(*) quest
      iusepbc=1
      if (lquest .gt. 0) call askyn(quest,lquest,1,-1,iusepbc,0,0)
      if (iusepbc .eq. 1) then
        call pbctype(ioppbc,npbc,ixyzhex,nonone)
        call pbcsize(ioppbc,edge,npbc)
        if (ioppbc .eq. 0) then
          call readimg(cell,ncell)
        else
          call crorgn(edge,edge_gen,ioppbc,3,ncell,cell,cellalt,
     -      ixyzhex,rinscr,rcirc)
        end if
        call prtcell(ioppbc,edge,edge_gen,0.0,vol,nw,1)
      end if
      if (ioppbc .eq. -1) iusepbc=0
      return
      end
