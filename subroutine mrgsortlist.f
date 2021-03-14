      subroutine mrgsortlist(list,it1,it2,it3,it4,it5,n)
      dimension list(n),it1(n),it2(n),it3(n),it4(n),it5(n)
c     Sort a list of positive integers
c     print *,'MRGS n=',n
c     print *,'LIST=',(list(i),i=1,n)
      call indexit(it1,1,n,0)
c     do i=1,n
c       t1(i)=list(i)
c     end do
      call mrgsrti(6,it1,list,n,it2,it3,it4,it5,n)
c     print *,'T=',(t1(i),i=1,n)
c     do i=1,n
c       list(i)=t1(i)
c     end do
c     print *,'LIST=',(list(i),i=1,n)
      return
      end
