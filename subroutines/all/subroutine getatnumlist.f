      subroutine getatnumlist(n,iatnum,ifchrg,ialist,icatlist,ixlist,
     -  nanos)
      dimension iatnum(n),ifchrg(n),ixlist(99),ialist(15),icatlist(15)
c     Get a list of atomic numbers that occured in this system
      call zeroiti(ixlist,0,99)
      call zeroiti(icatlist,0,15)
      nanos=0
      do i=1,n
        ifound=0
        do j=1,nanos
          if (iatnum(i) .eq. ialist(j)) then
            ifound=1
            go to 9031
          end if
        end do
9031    if (ifound .eq. 0) then
          if (nanos .lt. 15) then
            nanos=nanos+1
            ialist(nanos)=iatnum(i)
            ixlist(iatnum(i))=nanos
            icatlist(nanos)=ifchrg(i)
          else
            print *,'ERROR: Number of different elements exceeds',
     -        ' 15 - redimension the program'
            stop
          end if
        end if
      end do
c      write (6,7712) 'ialist',(ialist(i),i=1,nanos)
c      write (6,7712) 'icatlist',(icatlist(i),i=1,nanos)
c7712  format(1x,a,'=',15i3)
      return
      end
