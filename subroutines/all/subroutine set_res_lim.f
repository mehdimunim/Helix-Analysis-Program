      subroutine set_res_lim(ires,n,ifres,ilres,nres,index,iuseindex,
     -  maxres,maxn)
      dimension ires(maxn),ifres(maxres),ilres(maxres),index(maxn)
c     print *,'SET_RES_LIM n=',n
      if (iuseindex .eq. 0) call indexit(index,1,n,0)
      ifres(1)=1
      iresprev=ires(index(1))
      nres=1
      do ia=2,n
        if (ires(index(ia)) .ne. iresprev) then
          ilres(nres)=ia-1
          nres=nres+1
          ifres(nres)=ia
          iresprev=ires(index(ia))
        end if
      end do
      ilres(nres)=n
      return
      end
