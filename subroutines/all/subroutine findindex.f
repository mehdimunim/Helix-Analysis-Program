      subroutine findindex(lab,llab,q,lq,list,llist,nlist,ix,maxlablen,
     -  ihelp)
      character*(*) q,lab,list(nlist)
      dimension llist(nlist)
      ix=0
110   call getname(lab,llab,q(1:lq),lq,maxlablen,'',0,0,ihelp,0)
c     print *,'FINDINDEX llab=',llab
      if (llab .le. 1) return
      do i=1,nlist
        if (lab(1:llab) .eq. list(i)(1:llist(i))) ix=i
      end do
      if (ix .eq. 0) then
        print *,'No match found for ',lab(1:llab)
        go to 110
      end if
      return
      end
