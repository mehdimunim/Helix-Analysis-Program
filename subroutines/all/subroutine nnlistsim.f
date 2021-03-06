      subroutine nnlistsim(nfirst,n,c,nneig,ineig,indices,nbox,
     -  rcut,ifail,maxng,nnlistlen,maxbox,maxrec,LEVTEST)
      dimension nneig(n),ineig(maxng,nnlistlen),c(3,n),
     -  indices(maxbox,maxrec),nbox(maxrec)
      dimension nxyz(3),ix(3),xyzmin(3),xyzmax(3),centinp(3)
c     Set up neighbour list using linked cells
c     print *,'NNLS maxng,nnlistlen,maxbox=',maxng,nnlistlen,maxbox
      if (n .gt. nnlistlen) then
        print *,'Redimension the program with maxat > ',10*n
        ifail=1
        return
      end if
      div=rcut+0.1
      rcut2=rcut**2
      ifail=0
      if (n .lt. nfirst) return
      call extension(c,nneig,0,nfirst,n,xyzmin,xyzmax,centinp,0,0,v)
      call zeroiti(nneig,0,n)
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
     -    ' - increase maxrec to ',ngrid,' to avoid it'
        do k=1,3
          nxyz(k)=(xyzmax(k)-xyzmin(k))/div+1
        end do
        ngrid=nxyz(1)*nxyz(2)*nxyz(3)
      end do
      nx=nxyz(1)
      nxy=nxyz(1)*nxyz(2)
      if (LEVTEST .gt. 0)
     -  write (88,1000) div,nx,nxy,nxyz,xyzmin,xyzmax
      call zeroiti(nbox,0,ngrid)
      nboxmax=0
      do i=nfirst,n
         do k=1,3
           ix(k)=(c(k,i)-xyzmin(k))/div+1
         end do
c        Save i into the box represented by the indices ix(1-3)
         indexi=ix(1)+nx*(ix(2)-1)+nxy*(ix(3)-1)
         if (LEVTEST .gt. 1) write (88,1001) i,indexi,ix
         nbox(indexi)=nbox(indexi)+1
         if (nboxmax .lt. nbox(indexi)) nboxmax=nbox(indexi)
         if (nbox(indexi) .le. maxbox) then
           indices(nbox(indexi),indexi)=i
         else
           ifail=1
         end if
      end do
c     print *,'nboxmax=',nboxmax
      if (ifail .gt. 0) then
        print *,'Too many atoms in a box - ',
     -    'increase maxbox or MAXREC and recompile'
        return
      end if
c     Loop on the boxes
      do i1=1,nxyz(1)
        do i2=1,nxyz(2)
          do i3=1,nxyz(3)
            indexi=i1+nx*(i2-1)+nxy*(i3-1)
            ni=nbox(indexi)
            if (LEVTEST .gt. 2) write (88,1002) i1,i2,i3,indexi,ni
            if (ni .gt. 0) then
              do j1=max0(1,i1-1),i1
                do j2=max0(1,i2-1),min0(nxyz(2),i2+(i1-j1))
                  j3lim=i3
                  if (i1 .gt. j1 .or. i2. gt. j2) j3lim=j3lim+1
                  do j3=max0(1,i3-1),min0(j3lim,nxyz(3))
                    indexj=j1+nx*(j2-1)+nxy*(j3-1)
                    nj=nbox(indexj)
c                   ni, nj are the number of atoms in the box (i1,i2,i3) and
c                   (j1,j2,j3)
                    if (i1 .eq. j1 .and. i2 .eq. j2 .and.
     -                  i3 .eq. j3) then
                      ij=0
                    else
                      ij=1
                    end if
                    if (LEVTEST .gt. 2)
     -                write (88,1012) j1,j2,j3,indexj,nj
                    do jj=1,nj
                      j=indices(jj,indexj)
                      ii1=1
                      if (ij .eq. 0) ii1=jj+1
                      do ii=ii1,ni
                        i=indices(ii,indexi)
                        r2=dist2(c(1,i),c(1,j))
                        if (r2 .le. rcut2) then
c                         Bond found.
                          if (nneig(i) .lt. maxng .and.
     -                        nneig(j) .lt. maxng) then
                          nneig(i)=nneig(i)+1
                          ineig(nneig(i),i)=j
                          nneig(j)=nneig(j)+1
                          ineig(nneig(j),j)=i
                          else
                            print *,'Redimension the program for more',
     -                        ' neighbors in nnlistsim'
                            ifail=1
                            return
                          end if
                        end if
                      end do
                    end do
                  end do
                end do
              end do
            end if
          end do
        end do
      end do
      maxnn=0
      do i=nfirst,n
        if (nneig(i) .gt. maxnn) maxnn=nneig(i)
      end do
c     print *,'maxnn=',maxnn
      return
1000  format(' NNLISTSIM div=',f10.5,' nx,nxy=',2i4,
     -  ' nxyz=',3i4,/,' xyzmin=',3f10.5,' xyzmax=',3f10.5)
1001  format(' NNLISTSIM i,indexi=',2i6,' ix=',3i4)
1002  format(' NNLISTSIM i1,2,3=',3i4,' indexi=',i6,' ni=',i4)
1012  format(' NNLISTSIM j1,2,3=',3i4,' indexj=',i6,' nj=',i4)
      end
