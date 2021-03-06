      subroutine mrgsrti(iout,indexx,ivalue,n,ifa,ila,itemp1,itemp2,
     -  maxt)
c#    MMC routine 262 lstmod: 10/17/97
c*****Sort the indexx and ivalue array in the order of increasing ivalue(i)
      dimension indexx(n),ivalue(n)
      dimension ifa(maxt),ila(maxt),itemp1(maxt),itemp2(maxt)
      call mrglimtst(iout,n,maxt,ireturn)
      if (ireturn .eq. 1) return
c     Individual intervals consist of single elements
      nn=n
c     print *,'MRGST n,maxt=',n,maxt
      call indexit(ifa,1,n,0)
      call indexit(ila,1,n,0)
c     Merge pairs of intervals
11    nnpair=nn/2
      l=1
      do i=1,nnpair
        call mergelsti(indexx,ivalue,ifa(l),ila(l),ifa(l+1),ila(l+1),
     -    itemp1,itemp2,maxt)
        ifa(i)=ifa(l)
        ila(i)=ila(l+1)
        l=l+2
      end do
      if (2*nnpair .ne. nn) then
c       Take care of last (odd) interval
        call mergelsti(indexx,ivalue,ifa(nnpair),ila(nnpair),ifa(nn),
     -    ila(nn),itemp1,itemp2,maxt)
        ila(nnpair)=ila(nn)
      end if
      nn=nnpair
      if (nn .gt. 1) go to 11
      return
      end
