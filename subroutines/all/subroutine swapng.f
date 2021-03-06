      subroutine swapng(ineig,ibnd,ia,n1,n2,n,maxng)
      dimension ineig(maxng,n),ibnd(maxng,n)
c     Swap i-ia bond to the from n1 to n2 position
      if (n1 .lt. n2) then
        call swapi4(ineig(n1,ia),ineig(n2,ia))
c       ii=ineig(n1,ia)
c       ineig(n1,ia)=ineig(n2,ia)
c       ineig(n2,ia)=ii
        call swapi4(ibnd(n1,ia),ibnd(n2,ia))
c       ii=ibnd(n1,ia)
c       ibnd(n1,ia)=ibnd(n2,ia)
c       ibnd(n2,ia)=ii
      end if
      return
      end
