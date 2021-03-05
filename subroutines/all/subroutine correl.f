      subroutine correl(iout,res1,ic1,name1,lname1,av1,sd1,res2,ic2,
     -  name2,lname2,av2,sd2,corr,n,iavsdcalc,now6)
      dimension res1(2,n),res2(2,n)
      character*(*) name1,name2
      real*8 sum12,sum1,sum2,sum11,sum22
      if (iavsdcalc .gt. 0) then
c       Calculate and print avg & sd
        sum1=0.d0
        sum2=0.d0
        sum11=0.d0
        sum22=0.0d0
        do i=1,n
          sum1=sum1+res1(ic1,i)
          sum11=sum11+res1(ic1,i)**2
          sum2=sum2+res2(ic2,i)
          sum22=sum22+res2(ic2,i)**2
        end do
        av1=sum1/n
        av2=sum2/n
        sd1=dsqrt(dabs(sum11/n-av1**2))
        sd2=dsqrt(dabs(sum22/n-av2**2))
        write (6,1001) name1(1:lname1),av1,sd1
        write (6,1001) name2(1:lname2),av2,sd2
      end if
      sum12=0.d0
      do i=1,n
        sum12=sum12+res1(ic1,i)*res2(ic2,i)
      end do
      sum12=sum12/n
      corr=0.0
      if (sd1*sd2 .gt. 0.0) corr=(sum12-av1*av2)/(sd1*sd2)
      if (now6 .eq. 0)
     -  write (6,1000) name1(1:lname1),name2(1:lname2),corr
      write (iout,1000) name1(1:lname1),name2(1:lname2),corr
      return
1000  format(1x,a25,' - ',a25,' correlation=',f8.5)
1001  format(1x,a25,' Average=',f15.6,' S.D.=',f16.5)
      end
