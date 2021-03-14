      subroutine nnlistmpx(n,nanchorr,nanchorn,indexa,indexov,iatnum,
     -  nneig,nhbneig,c,isegno,indices,nbox,rmpxlim,it1,it2,it3,it4,it5,
     -  temp1,temp2,ifail,maxbox,maxrec)
      dimension indexa(n),indexov(n),iatnum(n),nneig(n),nhbneig(n),
     -  c(3,n),isegno(maxrec),indices(maxbox,maxrec),
     -  nbox(maxrec),it1(n),it2(n),it3(n),it4(n),it5(n),temp1(n),
     -  temp2(n)
      dimension nxyz(3),ix(3),xyzmin(3),xyzmax(3),centinp(3)
c     Set up mutually proximal pair list using linked cells
c     print *,'NNLISTMPX n=',n,' nanchorr,nanchorn=',nanchorr,nanchorn
c     print *,'NNLISTMPX n=',n,' maxbox,maxrec=',maxbox,maxrec
c     indexa(ia): Reference set atoms
c     indexov(ia): neighbour set atoms
      LEVTEST=0
      ifail=0
      call extension(c,nneig,0,1,n,xyzmin,xyzmax,centinp,0,1,v)
      div=rmpxlim
      if (rmpxlim .gt. 5.0) div=5.0
      rphmin=0.0
      call gridspace(c,isegno,1,n,xyzmin,xyzmax,div,
     -  nxyz,nx,nxy,nbox,ix,indices,ifail,1,LEVTEST,maxbox,maxrec)
      nboxes=nxyz(1)*nxyz(2)*nxyz(3)
      call zeroiti(it1,0,n)
c     Leave hydrogens out from the mpx search
      do ia=1,nanchorr
        if (iatnum(indexa(ia)) .ne. 1) it1(indexa(ia))=1
      end do
      do ia=1,nanchorn
        if (iatnum(indexov(ia)) .ne. 1) it1(indexov(ia))=2
        write (77,*) 'Ni=',ia,' Rix=',indexov(ia),
     -    ' IA=',iatnum(indexov(ia))
      end do
c     Sort the indices into contiguous list NONE/ANCHOR/NEIGHBOUR (0/1/2)
      do ib=1,nboxes
        nb=nbox(ib)
        call indexit(it2,1,nb,0)
        do ia=1,nb
          temp1(ia)=it1(indices(ia,ib))
        end do
        call mrgsrt(6,it2,temp1,nb,it3,it4,it5,temp2,n)
        call trnsfi(it5,indices(1,ib),nb)
        do ia=1,nb
          indices(ia,ib)=it5(it2(ia))
        end do
c       write (88,6941) 'it5',ib,(it5(ia),ia=1,nb)
c       write (88,6941) 'it1',ib,(it1(it5(ia)),ia=1,nb)
c       write (88,6941) 'ian',ib,(iatnum(it5(ia)),ia=1,nb)
c       write (88,6941) 'it2',ib,(it2(ia),ia=1,nb)
c       write (88,6941) 'ind',ib,(indices(ia,ib),ia=1,nb)
c6941   format(1x,a,i6,':',(20i6))
      end do
c     Find the limits of the NONE and ANHOR stretches
      do ib=1,nboxes
        if (nbox(ib) .gt. 0) then
          ia=1
          do while (it1(indices(ia,ib)) .eq. 0 .and. ia .lt. nbox(ib))
            ia=ia+1
          end do
          if (it1(indices(ia,ib)) .gt. 0) then
            if (it1(indices(ia,ib)) .eq. 1) then
              it3(ib)=ia
              do while (it1(indices(ia,ib)) .eq. 1 .and.
     -                  ia .lt. nbox(ib))
                ia=ia+1
              end do
              if (it1(indices(ia,ib)) .eq. 2) then
                it4(ib)=ia
              else
                it4(ib)=nbox(ib)+1
              end if
            else
c             No ones
              it3(ib)=nbox(ib)+1
              it4(ib)=ia
            end if
          else
c           All zeros
            it3(ib)=nbox(ib)+1
            it4(ib)=nbox(ib)+1
          end if
c         write (88,6791)ib,it3(ib),it4(ib),
c    -      (it1(indices(ia,ib)),ia=1,nbox(ib))
c         write (88,6792)ib,(indices(ia,ib),iatnum(indices(ia,ib)),
c    -      ia=1,nbox(ib))
c6791     format(i6,' it3,it4=',2i4,' it1(indices)=',(10i6))
c6792     format(i6,' indices,iatnos=',(10i6))
        end if
      end do
c     Anchor atoms: ia=it3(ib),it4(ib)-1
c     Neighbour atoms: ia=it4(ib),nbox(ib)
c     All heavy atoms are indexed in the grid
c     Loop on the boxes
      call zeroiti(nneig,0,n)
      call zeroiti(nhbneig,0,n)
      call zeroiti(it2,0,n)
c     it2(ia), temp2(ia): nearest neighbor atom and distance from anchor ia
c     and vice versa
      do ia=1,n
        temp2(ia)=999999.9
      end do
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
                    if (nj .gt. 0) then
                      if (indexi .eq. indexj) then
                        ij=0
                      else
                        ij=1
                      end if
                      if (LEVTEST .gt. 2)
     -                  write (88,1003) j1,j2,j3,indexj,nj
c                     Work with atom pairs in boxes indexi and indexj
                      do iaa=it3(indexi),it4(indexi)-1
                        ia=indices(iaa,indexi)
                        do jaa=it4(indexj),nj
                          ja=indices(jaa,indexj)
                          r2=dist2(c(1,ia),c(1,ja))
                          if (r2 .lt. temp2(ia)) then
                            temp2(ia)=r2
                            it2(ia)=ja
                          end if
                          if (r2 .lt. temp2(ja)) then
                            temp2(ja)=r2
                            it2(ja)=ia
                          end if
                        end do
                      end do
                      do iaa=it4(indexi),ni
                        ia=indices(iaa,indexi)
                        do jaa=it3(indexj),it4(indexj)-1
                          ja=indices(jaa,indexj)
                          r2=dist2(c(1,ia),c(1,ja))
                          if (r2 .lt. temp2(ia)) then
                            temp2(ia)=r2
                            it2(ia)=ja
                          end if
                          if (r2 .lt. temp2(ja)) then
                            temp2(ja)=r2
                            it2(ja)=ia
                          end if
                        end do
                      end do
                    end if
                  end do
                end do
              end do
            end if
          end do
        end do
      end do
      if (LEVTEST .gt. 0) then
        do ia=1,n
c9783     format(' ia=',i6,' it1,it2=',2i6,' r21,r22=',2f8.2)
          if (it1(ia) .eq. 1) then
            if (it2(ia) .gt. 0) then
              if (it2(it2(ia)) .eq. ia) then
                write (88,1004) ia,it2(ia),sqrt(temp2(ia))
              end if
            end if
          end if
        end do
      end if
      return
1002  format(' NNLISTMPX i1,2,3=',3i4,' indexi=',i6,' ni=',i4)
1003  format(' NNLISTMPX j1,2,3=',3i4,' indexj=',i6,' nj=',i4)
1004  format(' NNLISTMPX MPX pair:',2i6,' Dist=',f8.3,' A')
      end
