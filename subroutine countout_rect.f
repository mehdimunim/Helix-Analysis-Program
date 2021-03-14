      subroutine countout_rect(ifst,ilst,c,edge,nouttot,kmax,n)
      dimension c(3,n),edge(3)
      dimension nxyzp(3),nxyzm(3)
      call zeroiti(nxyzp,0,3)
      call zeroiti(nxyzm,0,3)
      do ia=ifst,ilst
        do k=1,3
          if (c(k,ia) .lt. 0.0) nxyzm(k)=nxyzm(k)+1
          if (c(k,ia) .gt. edge(k)) nxyzp(k)=nxyzp(k)+1
        end do
      end do
      nouttot=0
      noutmax=0
      kmax=0
      do k=1,3
        nouttot=nouttot+nxyzm(k)+nxyzp(k)
        if (nxyzp(k) .gt. noutmax) then
          noutmax=nxyzp(k)
          kmax=k
        end if
        if (nxyzm(k) .gt. noutmax) then
          noutmax=nxyzm(k)
          kmax=-k
        end if
      end do
c      write (6,1000) nxyzm,nxyzp
c      write (6,1001) ifst,ilst,nouttot,edge
c1000  format(' COUNTOUT nxyzm=',3i9,' nxyzp=',3i9)
c1001  format(' COUNTOUT ifst,ilst=',2i6,' nout=',i6,' e=',3f10.5)
      return
      end
