      subroutine getxyz(q1,lq1,q2,lq2,default,c,noneg,ihelp)
      character*(*) q1,q2
      dimension c(3)
      character*1 xyz(3)
      character*80 quiz
      data xyz /'X','Y','Z'/
      quiz(1:lq1)=q1(1:lq1)
      lquiz=lq1+1
      quiz(lquiz+1:lquiz+lq2)=q2(1:lq2)
      lquiz=lquiz+lq2
      do k=1,3
        quiz(lq1+1:lq1+1)=xyz(k)
        call getreal(quiz,lquiz,default,c(k),noneg,ihelp)
      end do
      return
      end
