      subroutine findshift(c,r,list,c0,del,is,nlist,rlambda,lmax,nzero,
     -  maxrec)
      dimension c(3,maxrec),r(maxrec),list(maxrec),c0(3),del(3)
c     Calculate the allowed shift of center, rlambda, along del
      nzero=0
      rl2=r(is)
      rl=sqrt(rl2)
      coslrl=0.0
      do k=1,3
        coslrl=coslrl+del(k)*(c(k,list(is))-c0(k))
      end do
c     cosl=cosl/rl
      rlambda=rl2
      do i=1,nlist
        ri2=r(i)
        ri=sqrt(ri2)
        cosiri=0.0
        do k=1,3
          cosiri=cosiri+del(k)*(c(k,list(i))-c0(k))
        end do
c       cosi=cosi/ri
c       if (rl*cosl .ne. ri*cosi) then
c         rlambdai=(rl2-ri2)/(2.0*(rl*cosl-ri*cosi))
        if (coslrl .ne. cosiri) then
          rlambdai=(rl2-ri2)/(2.0*(coslrl-cosiri))
        else
          rlambdai=0.0
        end if
c        write (6,7777) i,rl2,ri2,rl,ri,cosl,cosi,rlambdai
c7777    format(i3,' rli2=',2f8.3,' rli=',2f6.2,' cli=',2f8.5,
c     -    ' l=',f8.3)
        if (rlambdai .lt. rlambda .and. rlambdai .ge. 0.0) then
          rlambda=rlambdai
          lmax=i
        end if
        if (rlambdai .lt. 1.e-4 .and. rlambdai .ge. 0.0) then
          call swapl(list,r,nlist-nzero,lmax,-1,maxrec)
          lmax=nlist-nzero
          nzero=nzero+1
        end if
      end do
      return
      end
