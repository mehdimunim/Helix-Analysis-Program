      subroutine mrgsrt(iout,indexx,value,n,ifa,ila,itemp,temp,maxt)
c*****Sort the indexx and value array in the order of increasing value(i)
      dimension indexx(n),value(n)
      dimension ifa(maxt),ila(maxt),itemp(maxt),temp(maxt)
c     print *,'MRGSRT n,maxt=',n,maxt
      call mrglimtst(iout,n,maxt,ireturn)
      if (ireturn .eq. 1) return
      nn=n
      call indexit(ifa,1,n,0)
      call indexit(ila,1,n,0)
c     Merge pairs of intervals
11    nnpair=nn/2
      l=1
      do i=1,nnpair
        call mergelst(indexx,value,ifa(l),ila(l),ifa(l+1),ila(l+1),
     -    itemp,temp,maxt)
        ifa(i)=ifa(l)
        ila(i)=ila(l+1)
        l=l+2
      end do
      if (2*nnpair .ne. nn) then
c       Take care of last (odd) interval
        call mergelst(indexx,value,ifa(nnpair),ila(nnpair),ifa(nn),
     -    ila(nn),itemp,temp,maxt)
        ila(nnpair)=ila(nn)
      end if
      nn=nnpair
      if (nn .gt. 1) go to 11
      return
      end
