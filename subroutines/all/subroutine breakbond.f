      subroutine breakbond(i1,i2,n0,n,nneig,ineig,nneiga,nneigh,ian,
     -  ifail,maxng)
      dimension nneig(n),ineig(maxng,n),nneiga(n),nneigh(n),ian(n)
      ifail=0
      if (i1 .lt. n0 .or. i1 .gt. n .or.
     -    i2 .lt. n0 .or. i2 .gt. n) then
        print *,'ERROR: invalid atomindices: ',i1,i2
        ifail=1
        return
      end if
      idel=0
      do in=1,nneiga(i1)
        if (ineig(in,i1) .eq. i2) idel=in
      end do
      if (idel .eq. 0) then
        write (6,1000) i2,i1
        ifail=1
      else
        if (ian(ineig(idel,i1)) .eq. 1) nneigh(i1)=nneigh(i1)-1
        do in=idel+1,nneig(i1)
          ineig(in-1,i1)=ineig(in,i1)
        end do
        nneig(i1)=nneig(i1)-1
        nneiga(i1)=nneiga(i1)-1
      end if
      return
1000  format(' ERROR: atom ',i5,' is not on the neighbour list of atom',
     -  i6)
      end
