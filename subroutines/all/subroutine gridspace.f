      subroutine gridspace(c,isegno,nfirst,n,xyzmin,xyzmax,div,nxyz,
     -  nx,nxy,nbox,ix,indices,ifail,istop,LEVTEST,maxbox,maxrec)
      dimension c(3,n),isegno(n),nxyz(3),ix(3),xyzmin(3),xyzmax(3),
     -  nbox(maxrec),indices(maxbox,maxrec)
c     print *,'GRIDSPACE MAXBOX,MAXREC=',maxbox,maxrec,' DIV=',div
      do k=1,3
        nxyz(k)=(xyzmax(k)-xyzmin(k))/div+1
      end do
      nx=nxyz(1)
      nxy=nxyz(1)*nxyz(2)
      if (LEVTEST .gt. 0)
     -  write (88,1000) div,nx,nxy,nxyz,xyzmin,xyzmax
      ngrid=nxyz(1)*nxyz(2)*nxyz(3)
      do while (ngrid .gt. maxrec)
c       Increase gridsize to reduce the number of boxes under maxrec
        divo=div
        div=div*(float(ngrid)/float(maxrec))**(1.0/3.0) + 0.1
        print *,'Grid size increased from ',divo,' to ',div,
     -    ' - increase MAXREC to ',ngrid,' to avoid it'
        do k=1,3
          nxyz(k)=(xyzmax(k)-xyzmin(k))/div+1
        end do
        ngrid=nxyz(1)*nxyz(2)*nxyz(3)
        nx=nxyz(1)
        nxy=nxyz(1)*nxyz(2)
      end do
      if (LEVTEST .gt. 0)
     -  write (88,1000) div,nx,nxy,nxyz,xyzmin,xyzmax
      call zeroiti(nbox,0,ngrid)
      maxboxcount=0
c     print *,'NNLIST00 NFIRST,N=',nfirst,n
c     nboxmax=0
      do i=nfirst,n
        if (isegno(i) .ne. -1) then
          do k=1,3
            ix(k)=(c(k,i)-xyzmin(k))/div+1
          end do
c         Save i into the box represented by the indices ix(1-3)
          indexi=ix(1)+nx*(ix(2)-1)+nxy*(ix(3)-1)
          if (LEVTEST .gt. 1) write (88,1001) i,indexi,ix
          nbox(indexi)=nbox(indexi)+1
c         if (nboxmax .lt. nbox(indexi)) nboxmax=nbox(indexi)
          if (nbox(indexi) .le. maxbox) then
            indices(nbox(indexi),indexi)=i
          else
            print *,'Too many atoms in a box'
            if (istop .eq. 0) then
              print *,'Slow bond generator will be used for now'
              print *,'Increase MAXNEIG and recompile to avoid it'
              ifail=1
              return
            else
              print *,'Reduce distance threshold or'
              print *,'increase MAXNEIG and recompile'
              stop
            end if
          end if
        end if
      end do
c     print *,'NBOXMAX=',nboxmax
      return
1000  format(' GRIDSPACE div=',f7.2,' nx,nxy=',2i8,
     -  ' nxyz=',3i4,/,' xyzmin=',3f10.5,' xyzmax=',3f10.5)
1001  format(' GRIDSPACE i,indexi=',2i6,' ix=',3i4)
      end
