      subroutine breakbond0(i1,i2,n,nneig,ineig,ifail,maxng)
      dimension nneig(n),ineig(maxng,n)
      ifail=0
c     print *,'BREAKBIND0 nneig(i1,i2),n=',nneig(i1),nneig(i2),n
      if (i1 .lt. 1 .or. i1 .gt. n .or.
     -    i2 .lt. 1 .or. i2 .gt. n) then
        print *,'ERROR: invalid atomindices: ',i1,i2
        ifail=1
        return
      end if
      idel=0
      do in=1,nneig(i1)
        if (ineig(in,i1) .eq. i2) idel=in
      end do
      if (idel .eq. 0) then
        write (6,1000) i2,i1
        ifail=1
      else
        do in=idel+1,nneig(i1)
          ineig(in-1,i1)=ineig(in,i1)
        end do
        nneig(i1)=nneig(i1)-1
      end if
      idel=0
      do in=1,nneig(i2)
        if (ineig(in,i2) .eq. i1) idel=in
      end do
      if (idel .eq. 0) then
        write (6,1000) i1,i2
        ifail=1
      else
        do in=idel+1,nneig(i2)
          ineig(in-1,i2)=ineig(in,i2)
        end do
        nneig(i2)=nneig(i2)-1
      end if
      return
1000  format(' ERROR: atom ',i5,' is not on the neighbour list of atom',
     -  i6)
      end
