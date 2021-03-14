      subroutine bondlenstat(c,n,iatnum,nneig,ineig,iatnm2,nconfig,
     -  nanos,ixlist,ialist,iwbls,inpcrdtyp,ioins,ic1,ic2,ir1,ir2,
     -  line,index,radtodeg,maxrepconf,maxng,maxrec)
      dimension c(3,n),iatnum(n),nneig(n),ineig(maxng,n),index(n),
     -  ixlist(99),ialist(15)
      character* 132 line(maxrec)
      character*2 iatnm2(99)
      character*1 sp
      dimension dmin(15,15),dav(15,15),dsd(15,15),ijmin(2,15,15),
     -  nij(15,15),dmax(15,15),ijmax(2,15,15),dd(15),anmin(15,15,15),
     -  anav(15,15,15),ansd(15,15,15),ijkmin(3,15,15,15),nijk(15,15,15),
     -  anmax(15,15,15),ijkmax(3,15,15,15),anan(15,15)
      real*8 cosa
      data ixl /0/,jxl /0/,sp /' '/
c     print *,'BLS ic1,ic2,ir1,ir2=',ic1,ic2,ir1,ir2
c     Bond length statistics first
      do i=1,nanos
        do j=1,i
c         j<=i
          nij(i,j)=0
          dav(i,j)=0.0
          dsd(i,j)=0.0
          dmax(i,j)=0.0
          dmin(i,j)=10000.0
          ijmin(1,i,j)=0
          ijmax(1,i,j)=0
          ijmin(2,i,j)=0
          ijmax(2,i,j)=0
        end do
      end do
      do i=1,n
        ixll=ixlist(iatnum(i))
        do jj=1,nneig(i)
          j=ineig(jj,i)
          jxll=ixlist(iatnum(j))
          dd(jj)=sqrt((c(1,i)-c(1,j))**2+(c(2,i)-c(2,j))**2+
     -       (c(3,i)-c(3,j))**2)
          if (jxll .le. ixll) then
            ixl=ixll
            jxl=jxll
          else
            ixl=jxll
            jxl=ixll
          end if
          dav(ixl,jxl)=dav(ixl,jxl)+dd(jj)
          dsd(ixl,jxl)=dsd(ixl,jxl)+dd(jj)**2
          nij(ixl,jxl)=nij(ixl,jxl)+1
          if (dmax(ixl,jxl) .lt. dd(jj)) then
            dmax(ixl,jxl)=dd(jj)
            ijmax(1,ixl,jxl)=i
            ijmax(2,ixl,jxl)=j
          end if
          if (dmin(ixl,jxl) .gt. dd(jj)) then
            dmin(ixl,jxl)=dd(jj)
            ijmin(1,ixl,jxl)=i
            ijmin(2,ixl,jxl)=j
          end if
        end do
        if (inpcrdtyp .le. ioins) then
          write (iwbls,2046) i,line(index(i))(ic1:ic2),
     -      line(index(i))(ir1:ir2),
     -      (sp,ineig(j,i),line(index(ineig(j,i)))(ic1:ic2),
     -      line(index(ineig(j,i)))(ir1:ir2),dd(j),j=1,nneig(i))
        else
          write (iwbls,2048) i,iatnm2(iatnum(i)),
     -      (sp,ineig(j,i),iatnm2(iatnum(ineig(j,i))),
     -      dd(j),j=1,nneig(i))
        end if
      end do
      if (nconfig .le. maxrepconf) write (6,2057)
      write (iwbls,2057)
      do i=1,nanos
        do j=1,i
          if (nij(i,j) .gt. 0) dav(i,j)=dav(i,j)/nij(i,j)
          dsd(i,j)=sqrt(abs(dsd(i,j)/nij(i,j)-dav(i,j)**2))
          if (dmin(i,j) .eq. 10000.0) dmin(i,j)=0.0
        end do
        do j=1,i
          if (ijmin(1,i,j) .gt. 0) then
            if (nconfig .le. maxrepconf)
     -       write (6,2047) iatnm2(ialist(i)),
     -         iatnm2(ialist(j)),dmin(i,j),ijmin(1,i,j),ijmin(2,i,j),
     -         dav(i,j),dsd(i,j),dmax(i,j),ijmax(1,i,j),ijmax(2,i,j)
            write (iwbls,2047) iatnm2(ialist(i)),
     -       iatnm2(ialist(j)),dmin(i,j),ijmin(1,i,j),ijmin(2,i,j),
     -       dav(i,j),dsd(i,j),dmax(i,j),ijmax(1,i,j),ijmax(2,i,j)
          end if
        end do
      end do
c     Bond angle statistics next
      do i=1,nanos
        do j=1,nanos
          do k=1,j
c           k<=j
            nijk(i,j,k)=0
            anav(i,j,k)=0.0
            ansd(i,j,k)=0.0
            anmax(i,j,k)=0.0
            anmin(i,j,k)=10000.0
            ijkmin(1,i,j,k)=0
            ijkmax(1,i,j,k)=0
            ijkmin(2,i,j,k)=0
            ijkmax(2,i,j,k)=0
            ijkmin(3,i,j,k)=0
            ijkmax(3,i,j,k)=0
          end do
        end do
      end do
      do i=1,n
        ixll=ixlist(iatnum(i))
        do jj=1,nneig(i)
          do kk=1,jj-1
            j=ineig(jj,i)
            k=ineig(kk,i)
            jxll=ixlist(iatnum(j))
            kxll=ixlist(iatnum(k))
            ixl=max0(ixll,jxll,kxll)
            kxl=min0(ixll,jxll,kxll)
            if (ixll .ne. ixl .and. ixll .ne. kxl) jxl=ixll
            if (jxll .ne. ixl .and. jxll .ne. kxl) jxl=jxll
            if (kxll .ne. ixl .and. kxll .ne. kxl) jxl=kxll
            rj=0.0
            rk=0.0
            rjk=0.0
            do l=1,3
              dj=c(l,i)-c(l,j)
              dk=c(l,i)-c(l,k)
              rj=rj+dj*dj
              rk=rk+dk*dk
              rjk=rjk+dj*dk
            end do
            cosa=dble(rjk)/sqrt(rj*rk)
            anan(jj,kk)=radtodeg*dacoscheck(cosa,ccc,1,6,'BONDLENSTAT')
            anav(ixl,jxl,kxl)=anav(ixl,kxl,jxl)+anan(jj,kk)
            ansd(ixl,jxl,kxl)=ansd(ixl,kxl,jxl)+anan(jj,kk)**2
            nijk(ixl,jxl,kxl)=nijk(ixl,kxl,jxl)+1
            if (anmax(ixl,jxl,kxl) .lt. anan(jj,kk)) then
              anmax(ixl,jxl,kxl)=anan(jj,kk)
              ijkmax(1,ixl,jxl,kxl)=i
              ijkmax(2,ixl,jxl,kxl)=j
              ijkmax(3,ixl,jxl,kxl)=k
            end if
            if (anmin(ixl,jxl,kxl) .gt. anan(jj,kk)) then
              anmin(ixl,jxl,kxl)=anan(jj,kk)
              ijkmin(1,ixl,jxl,kxl)=i
              ijkmin(2,ixl,jxl,kxl)=j
              ijkmin(3,ixl,jxl,kxl)=k
            end if
            if (inpcrdtyp .le. ioins) then
              write (iwbls,2050) j,line(index(j))(ic1:ic2),
     -          line(index(j))(ir1:ir2),i,line(index(i))(ic1:ic2),
     -          line(index(i))(ir1:ir2),k,line(index(k))(ic1:ic2),
     -          line(index(k))(ir1:ir2),anan(jj,kk)
            else
              write (iwbls,2051) j,iatnm2(iatnum(j)),
     -          i,iatnm2(iatnum(i)),k,iatnm2(iatnum(k)),anan(jj,kk)
            end if
          end do
        end do
      end do
      if (nconfig .le. maxrepconf) write (6,2058)
      write (iwbls,2058)
      do i=1,nanos
        do j=1,nanos
          do k=1,j
            if (nijk(i,j,k) .gt. 0) then
              anav(i,j,k)=anav(i,j,k)/nijk(i,j,k)
              ansd(i,j,k)=sqrt(abs(ansd(i,j,k)/nijk(i,j,k)-
     -          anav(i,j,k)**2))
            end if
            if (anmin(i,j,k) .eq. 10000.0) anmin(i,j,k)=0.0
            if (ijkmin(1,i,j,k) .gt. 0) then
              if (nconfig .le. maxrepconf)
     -          write (6,2049) iatnm2(ialist(j)),iatnm2(ialist(i)),
     -            iatnm2(ialist(k)),anav(i,j,k),ansd(i,j,k),
     -            ijkmin(2,i,j,k),ijkmin(1,i,j,k),ijkmin(3,i,j,k),
     -            anmin(i,j,k),ijkmax(2,i,j,k),ijkmax(1,i,j,k),
     -            ijkmax(3,i,j,k),anmax(i,j,k)
              write (iwbls,2049) iatnm2(ialist(j)),iatnm2(ialist(i)),
     -          iatnm2(ialist(k)),anav(i,j,k),ansd(i,j,k),
     -          ijkmin(2,i,j,k),ijkmin(1,i,j,k),ijkmin(3,i,j,k),
     -          anmin(i,j,k),ijkmax(2,i,j,k),ijkmax(1,i,j,k),
     -          ijkmax(3,i,j,k),anmax(i,j,k)
            end if
          end do
        end do
      end do
      return
2046  format(i5,' From ',a4,1x,a4,a1,'To ',i5,1x,a4,1x,a4,f6.3,';',a1,
     -                            ' To ',i5,1x,a4,1x,a4,f6.3,';',a1,
     -  /,(20x,2(' To ',i5,1x,a4,1x,a4,f6.3,';',a1)))
2047  format(1x,a2,' - ',a2,f6.3,' (',i5,'-',i5,') < ',f6.3,' (av) +/-',
     -  f4.3,' (SD) < ',f6.3,' (',i5,'-',i5,')')
2048  format(i5,' From ',a2,1x,' To ',i5,1x,a2,1x,f6.3,';',a1,
     -  ' To ',i5,1x,a2,1x,f6.3,';',a1,' To ',i5,1x,a2,1x,f6.3,';',a1,
     -  /,(15x,3(' To ',i5,1x,a2,1x,f6.3,';',a1)))
2049  format(' Angle ',a2,'-',a2,'-',a2,' Average=',f8.2,' SD=',f8.2,
     -  ' deg',/,6x,'Min (',i5,'-',i5,'-',i5,')=',f8.2,
     -  ' Max (',i5,'-',i5,'-',i5,')=',f8.2)
2050  format(' Angle ',i5,' (',a4,1x,a4,') - ',i5,' (',a4,1x,a4,') - ',
     -  i5,' (',a4,1x,a4,')=',f10.5,' deg')
2051  format(' Angle ',i5,' (',a2,') - ',i5,' (',a2,') - ',
     -  i5,' (',a2,')=',f10.5,' deg')
2057  format(' Minimum, average and maximum of bond-type lengths')
2058  format(' Minimum, average and maximum of angle-type values')
      end
