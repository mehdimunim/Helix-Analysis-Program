      subroutine pickname(query,lquery,list,llist,nlist,ix)
      character*(*) query,list(nlist)
      dimension llist(nlist)
      character*80 line
c     print *,'PICKNAME nlist=',nlist
      maxlen=0
      do i=1,nlist
        if (llist(i) .gt. maxlen) maxlen=llist(i)
      end do
      if (maxlen .gt. 76) then
        print *,'PROGRAM ERROR: list item is longer than 76 chars'
        stop
      end if
      do i=1,nlist
        line(1:llist(i))=list(i)(1:llist(i))
        call blankout(line,llist(i)+1,maxlen+1)
        write (line(maxlen+2:maxlen+3),1001) i
        write (6,1000) line(1:maxlen+3)
      end do
      call getint(query,lquery,0,1,nlist,ix,000)
      return
1000  format(1x,a)
1001  format(i2)
      end
