      subroutine nnlist00(nfirst,n,nslt,islvw,iatnum,ifchrg,c,nneig,
     -  nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,npneig,line,
     -  irescol1,irescol2,inamcol1,inamcol2,index,maxng,hblimfac,angmin,
     -  indices,nbox,ixres,isegno,ifail,nframe,radtodeg,maxbox,maxrec,
     -  LEVTEST)
      dimension nneig(n),ineig(maxng,n),iatnum(n),ifchrg(n),c(3,n),
     -  nhbneig(n),nhneig(n),nnneig(n),ncneig(n),nsneig(n),
     -  npneig(n),index(n),indices(maxbox,maxrec),nbox(maxrec),ixres(n),
     -  isegno(n)
      common /connatdat/ ramax(99),ramax2(99),hlimfac,ianfg(99),
     -  namfcg(100),nrmw
      character* 132 line(maxrec)
      character*4 atnami,atnamj
      character*8 resnami,resnamj
      dimension nxyz(3),ix(3),xyzmin(3),xyzmax(3),centinp(3),
     -  nh12(4)
c     Set up neighbour list using linked cells
c     nhneig, ncneig,nnneig,nsneig,npneig:
c     Number of H, C, N, S, P neighbours, resp.
c     print *,'NNLIST00 nfirst,n,nslt,islvw=',nfirst,n,nslt,islvw
      ifail=0
      if (n .lt. nfirst) return
c     Find minimum and maximum atomic radius and extents of the molecule
      ra2max=ramax2(iatnum(1))
      ra2min=ra2max
      do i=nfirst,n
        if (ra2min .gt. ramax2(iatnum(i))) ra2min=ramax2(iatnum(i))
        if (ra2max .lt. ramax2(iatnum(i))) ra2max=ramax2(iatnum(i))
c       write (77,9877) i,line(index(i))(inamcol1:inamcol2),
c    -    iatnum(i),ramax2(iatnum(i)),ra2max
c9877   format(' NNLIST0',i6,a,' ia',i3,' ramax2=',f5.2,' ra2max=',f5.2)
      end do
      ramn=sqrt(ra2min)
      ramx=sqrt(ra2max)
      if (ramx .lt. 0.01) then
        if (nrmw .eq. 0)
     -    print *,'WARNING: no atomic number could be deduced'
        ramx=0.5
        nrmw=nrmw+1
      end if
      div=ramx*hblimfac+0.01
c     print *,'DIV=',div,' RAMX=',ramx
      call extension(c,nneig,0,nfirst,n,xyzmin,xyzmax,centinp,0,0,v)
      call gridspace(c,isegno,nfirst,n,xyzmin,xyzmax,div,
     -  nxyz,nx,nxy,nbox,ix,indices,ifail,0,LEVTEST,maxbox,maxrec)
      nmb=0
      if (ifail .gt. 0) return
c     Loop on the boxes
      do i1=1,nxyz(1)
        do i2=1,nxyz(2)
          do i3=1,nxyz(3)
            indexi=i1+nx*(i2-1)+nxy*(i3-1)
            ni=nbox(indexi)
            if (nmb .lt. ni) nmb=ni
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
     -                write (88,1012) j1,j2,j3,indexj,nj
                    do jj=1,nj
                      j=indices(jj,indexj)
                      if (inamcol2 .ge. inamcol1) then
                        resnamj='     '
                        resnamj=line(index(j))(irescol1:irescol2)
                        atnamj=line(index(j))(inamcol1:inamcol2)
                      end if
                      ii1=1
                      if (ij .eq. 0) ii1=jj+1
                      do ii=ii1,ni
                        i=indices(ii,indexi)
                        numhyd=0
                        if (iatnum(j) .eq. 1) numhyd=1
                        if (iatnum(i) .eq. 1) numhyd=numhyd+1
c                       write (78,*) 'NNLIST00 i,j=',i,j,' numhyd=',numhyd
                        if (numhyd .lt. 2 .and. ramax2(iatnum(i))*
     -                      ramax2(iatnum(j)) .gt. 0.0) then
c                         At most one hydrogen
                          r2=dist2(c(1,i),c(1,j))
                          rlm=amax1(ramax2(iatnum(i)),ramax2(iatnum(j)))
                          if (numhyd .eq. 1) rlm=rlm*hlimfac
                          if (r2 .le. rlm .and.
     -                        ifchrg(i)+ifchrg(j) .eq. 0) then
c                           Bond found. Bonds to electron, charge or lone pair
c                           are not saved in the 'real' atoms' list here.
                            if (i .le. nslt .or. j .le. nslt) then
c                             Solvent-solvent connectivity is already done
                              call savebond(i,j,iatnum,nneig,ineig,
     -                          nhbneig,nhneig,ncneig,nnneig,npneig,
     -                          nsneig,resnamj,atnamj,maxng,n,0,
     -                          inamcol1,inamcol2,LEVTEST)
                              resnami='     '
                              resnami=line(index(i))(irescol1:irescol2)
                              atnami=line(index(i))(inamcol1:inamcol2)
                              call savebond(j,i,iatnum,nneig,ineig,
     -                          nhbneig,nhneig,ncneig,nnneig,npneig,
     -                          nsneig,resnami,atnami,maxng,n,0,
     -                        inamcol1,inamcol2,LEVTEST)
                            end if
                          else if (numhyd .eq. 1) then
c                           No bond, but may be hydrogen bond
                            if (iatnum(i)*iatnum(j) .ne. 6 .and.
     -                          ifchrg(i)+ifchrg(j) .eq. 0) then
c                             No H-C or H-+
                              rlimhb2=(rlm/hlimfac)*hblimfac**2
c                             slv-slv bond only allowed for water bridge calc
c??                           Hbond with C between slt-slv?? - removed!!
                              if (islvw .eq. 2 .or.
     -                            i .le. nslt .or. j .le. nslt)
     -                          call maybehbond(r2,i,j,nneig,nhbneig,
     -                            ineig,inamcol1,inamcol2,rlimhb2,
     -                            atnami,resnami,n,maxng)
c                               if (i .eq.  900 .and. j .eq. 3142 .or.
c    -                              j .eq.  900 .and. i .eq. 3142)    
c    -write (6,8976) ,i,j,r2,rlimhb2,ramax2(iatnum(i)),ramax2(iatnum(j))
c8976 format('i,j=',2i6,' R2=',f6.2,'RLIMHB2=',f6.2,' ramax2:',2f7.2) 
                            end if
                          else if (numhyd .eq. 0 .and. islvw .gt. 0)then
                            if (i .gt. nslt .and. ifchrg(j) .gt. 0 .or.
     -                          j .gt. nslt .and. ifchrg(i) .gt. 0) then
                              rlimhb2=rlm*hblimfac**2
                              call maybehbond(r2,i,j,nneig,nhbneig,
     -                          ineig,inamcol1,inamcol2,rlimhb2,atnami,
     -                          resnami,n,maxng)
c                        print *,'MAYBE H-+:',i,j,nhbneig(i),nhbneig(j)
c                        print *,'MAYBE H-+: r2=',r2,' rlimhb2=',rlimhb2
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
      if (LEVTEST .gt. 0)
     -  write (6,2001) nmb,100.0*float(nmb)/float(maxbox)
c     Now screen the H-bonds for the angle and for C-H...X bond
      do ia=nfirst,n
        if (iatnum(ia) .eq. 1) then
c         Atom ia is always the H of the H bond
c         ihb0 is the heavy atom of the donor H
          call get_heavyat(ia,nneig,ineig,ixres,nframe,ihb0,maxng,
     -      maxrec)
          nhbng=nhbneig(ia)
          do ja=1,nhbneig(ia)
            ihb=ineig(maxng+1-ja,ia)
            if (ihb .gt. 0) then
              ihbdel=0
              if (ihb0 .gt. 0) call checkhbclose(c,n,ia,ihb,
     -          'hydrogen bonded',15,ihbdel)
c             ihb0 is the heavy atom of the donor H
c             ihb0=0
c             nng=nneig(ia)
c             do while (nng .gt. 0)
c               ihb0=ineig(nng,ia)
c               if (iatnum(ineig(nng,ia)) .eq. 1)  then
c                 ihb0=ineig(nng,ia)
c                 nng=0
c               end if
c               nng=nng-1
c             end do
              if (ihb0 .gt. 0 .and. ihbdel .eq. 0) then
                call checkhbclose(c,n,ia,ihb0,
     -            'hydrogen bonded',15,ihbdel)
                call checkhbclose(c,n,ihb0,ihb,
     -            'hydrogen-bond donor bond',24,ihbdel)
              end if
              if (ihb0 .gt. 0 .and. ihbdel .eq. 0) then
                call angdistw(c(1,ia),c(1,ihb),c(1,ihb0),rHB,rb,rab,ang)
                ang=ang*radtodeg
              else
c               Hydrogen had no heavy atom bonded to it - discard
                ang=angmin-1.0
              end if
              if (ang .lt. angmin .or. iatnum(ihb0) .eq. 6) then
c               Remove ihb from the HB list of ia
c               write (77,*) 'rHB,rb,rab=',rHB,rb,rab
                ineig(maxng+1-ja,ia)=0
c               Remove ia from the HB list of ihb
                nremia=0
                do ii=maxng+1-nhbneig(ihb),maxng
                  if (ineig(ii,ihb) .eq. ia) then
                    ineig(ii,ihb)=0
                    nremia=1
                  end if
                end do
                if (nremia .eq. 0) then
                  print *,'PROGRAM ERROR: atom ',ia,
     -            ' is not on the HB list of atom',ihb
                  write (6,*) '(trying to remove) ',ihb,' ihb:',
     -              (ineig(maxng-in+1,ihb),in=1,nhbneig(ihb))
                end if
              end if
            end if
          end do
        else if (ifchrg(ia) .gt. 0 .and. islvw .gt. 0) then
c         Cation - only H-bonds to water oxygen
          nhbng=nhbneig(ia)
          do ja=1,nhbneig(ia)
            ihb=ineig(maxng+1-ja,ia)
            if (ihb .gt. 0) then
c             Not dropped yet
              idrop=0
              if (iatnum(ihb) .eq. 1) then
                idrop=1
              else
                if (ihb .le. nslt) print *,'PROGRAM ERROR: atom ',ia,
     -            ' is a cation but is ','H-bonded to a non-solvent'
                if (iatnum(ihb) .ne. 8) print *,'PROGRAM ERROR: atom ',
     -            ia,' is H-bonded to a cation but is not an oxygen'
                nnh=0
                do in=1,nneig(ihb)
                  jn=ineig(in,ihb)
                  if (iatnum(jn) .eq. 1) then
                    nnh=nnh+1
                    nh12(nnh)=jn
                  end if
                end do
                if (nnh .ne. 2) then
                  if (nframe .eq. 0) write (6,2000) ihb,nnh
                  if (nframe .gt. 0) write (6,2000) ihb,nnh,' ',nframe
                end if
c               write (77,*) 'ia,ihb,nh12=',ia,ihb,nh12
                call angdistw(c(1,ihb),c(1,ia),c(1,nh12(1)),rHB,rb,
     -            rab,ang)
                ang=ang*radtodeg
                if (ang .lt. angmin) then
                  call angdistw(c(1,ihb),c(1,ia),c(1,nh12(2)),rHB,rb,
     -              rab,ang)
                  ang=ang*radtodeg
                end if
                if (ang .lt. angmin) idrop=2
              end if
              if (idrop .gt. 0) then
C               Drop
                print *,'DROP cation bond ',ia,' - ',ihb,
     -            ' idrop=',idrop,' ang=',ang
                ineig(maxng+1-ja,ia)=0
c               Remove ia from the HB list of ihb
                nremia=0
                do ii=maxng+1-nhbneig(ihb),maxng
                  if (ineig(ii,ihb) .eq. ia) then
                    ineig(ii,ihb)=0
                    nremia=1
                  end if
                end do
                if (nremia .eq. 0) then
                  print *,'PROGRAM ERROR: atom ',ia,
     -            ' is not on the HB list of atom',ihb
                  write (6,*) '(trying to remove) ',ihb,' ihb:',
     -              (ineig(maxng-in+1,ihb),in=1,nhbneig(ihb))
                end if
              end if
            end if
          end do
        end if
      end do
c     Now condense the list to fill in the zeros
      do ia=nfirst,n
        ndel=0
        do ja=1,nhbneig(ia)
          if (ineig(maxng-ja+1,ia) .eq. 0) then
            ndel=ndel+1
          else
            ineig(maxng-ja+1+ndel,ia)=ineig(maxng-ja+1,ia)
          end if
        end do
        nhbneig(ia)=nhbneig(ia)-ndel
      end do
      call checkhblist(n,ineig,nhbneig,maxng)
      return
1002  format(' NNLIST00 i1,2,3=',3i4,' indexi=',i6,' ni=',i4)
1012  format(' NNLIST00 j1,2,3=',3i4,' indexj=',i6,' nj=',i4)
2000  format(' WARNING: water oxygen ',i6,' has',i3,' hydrogens bonded',
     -  ' to it',a,' Nframe=',i6)
2001  format(' Geometry-based linked-cell NN search, max # of ',
     -  'atoms/cell=',i3,' (',f5.1,' %)')
      end
