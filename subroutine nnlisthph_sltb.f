      subroutine nnlisthph_sltb(n,ianchor2,iselfanc,indexa,nneig,
     -  n14neig,nhbneig,ineig,c,rhphmax,isltb,isegno,nosameseg,bondlab,
     -  lbondlab,indices,nbox,ifail,maxng,maxbox,maxrec)
      dimension indexa(n),nneig(n),n14neig(n),nhbneig(n),
     -  ineig(maxng,n),c(3,n),isegno(n),indices(maxbox,maxrec),
     -  nbox(maxrec)
      character*(*) bondlab
      dimension nxyz(3),ix(3),xyzmin(3),xyzmax(3),centinp(3)
c     Set up hydrophobic bond list using linked cells
c     print *,'NNLISTHPH_SLTB n,nosameseg,iselfanc,rhphmax=',
c    -  n,nosameseg,iselfanc,rhphmax
c     indexa(ia)=1 or 2: Anchor atom
c     indexa(ia)=-1 or -2: Non-anchor atom, but can form HPH bond/salt bridge
c     iabs(indexa(ia))=1: +; =2: -
      LEVTEST=0
      rhphmax2=rhphmax**2
      ifail=0
      call extension(c,nneig,0,1,n,xyzmin,xyzmax,centinp,0,1,v)
      div=amin1(5.001,rhphmax+0.001)
      rphmin=0.0
      call gridspace(c,isegno,1,n,xyzmin,xyzmax,div,
     -  nxyz,nx,nxy,nbox,ix,indices,ifail,1,LEVTEST,maxbox,maxrec)
c     All hydrophobic carbons are indexed in the grid
c     Loop on the boxes
      call zeroiti(nneig,0,n)
      call zeroiti(nhbneig,0,n)
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
                    if (indexi .eq. indexj) then
                      ij=0
                    else
                      ij=1
                    end if
                    if (LEVTEST .gt. 2)
     -                write (88,1003) j1,j2,j3,indexj,nj
                    do jj=1,nj
                      ja=indices(jj,indexj)
                      ii1=1
                      if (ij .eq. 0) ii1=jj+1
                      do ii=ii1,ni
                        ia=indices(ii,indexi)
                        idoit=0
                        if (indexa(ia) .gt. 0 .and.
     -                      indexa(ja) .gt. 0) then
                          if (iselfanc .eq. 1) idoit=1
                        else if (indexa(ia)*indexa(ja) .lt. 0) then
                          if (ianchor2 .eq. 0) idoit=1
                        end if
                        if (nosameseg .eq. 1 .and.
     -                      isegno(ia) .eq. isegno(ja)) idoit=0
                        if (isltb .eq. 1) then
c                         Salt bridge - exclude carbons
                          if (idoit .eq. 1) then
                            if (iabs(indexa(ia)) .eq. iabs(indexa(ja)))
     -                        idoit=0
                          end if
                        end if
                        if (idoit .eq. 1) then
c                         Make sure ia and ja far enough apart in topology
                          do in=1,n14neig(ia)
                            if (ja .eq. ineig(in,ia)) idoit=0
                          end do
c                         write (77,7211) ia,ja,idoit,
c    -                      (ineig(in,ia),in=1,n14neig(ia))
c7211                      format(' ia,ja,idoit=',3i5,' in(ia)=',6i5)
                        end if
                        if (idoit .eq. 1) then
                          r2=dist2(c(1,ia),c(1,ja))
                          if (r2 .le. rhphmax2) then
c                           Bond found
                            if (nhbneig(ia) .eq. maxng .or.
     -                          nhbneig(ia) .eq. maxng) then
                              write (6,2000) bondlab(1:lbondlab),maxng
                              ifail=1
                              return
                            else
                              nhbneig(ia)=nhbneig(ia)+1
                              ineig(maxng-nhbneig(ia)+1,ia)=ja
                              nhbneig(ja)=nhbneig(ja)+1
                              ineig(maxng-nhbneig(ja)+1,ja)=ia
                            end if
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
      if (LEVTEST .gt. 0) then
        do ia=1,n
          if (nhbneig(ia) .gt. 0) then
            write (88,1004) ia,(ineig(maxng-ja+1,ia),ja=1,nhbneig(ia))
          end if
        end do
      end if
      return
1002  format(' NNLISTHPH i1,2,3=',3i4,' indexi=',i6,' ni=',i4)
1003  format(' NNLISTHPH j1,2,3=',3i4,' indexj=',i6,' nj=',i4)
1004  format(' NNLISTHPH ia=',i4,' nhblist:',100i4)
2000  format(' ERROR: Maximum number of ',a,' neighbors (',i3,') is ',
     -  'exceeded')
      end
