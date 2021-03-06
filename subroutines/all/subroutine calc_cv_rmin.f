      subroutine calc_cv_rmin(c,ian,nslt,nslt1,nslt2,n,numsolv,naslv,
     -  iarepslv,spacing,rsltmax,rcvmax,cvmin,rnearsq,ianear,cv,cvsums,
     -  ncvsums,nconfig,indexs,ifirst,ilast,itemp1,itemp2,itemp3,iout,
     -  maxrec)
      dimension c(3,n),ian(n),rnearsq(maxrec),ianear(maxrec),
     -  cv(maxrec),ncvsums(maxrec),indexs(maxrec),ifirst(maxrec),
     -  ilast(maxrec),itemp1(maxrec),itemp2(maxrec),itemp3(maxrec)
      real*8 cvsums(3,maxrec)
      dimension edge(3),corner(3),ecell(3),xyzmin(3),xyzmax(3),nxyz(3),
     -  ixyz(3),ixyz1(3),ixyz2(3),d(3),cslv(3)
      real*8 cvsum(3)
c     print *,'CALC_CV_ numsolv,naslv,iarepslv=',numsolv,naslv,iarepslv
c     print *,'RSLTMAX,RCVMAX,CVMIN=',rsltmax,rcvmax,cvmin
c     print *,'CALC_CV nslt1,nslt2,nslt=',nslt1,nslt2,nslt
      call cellpart(c,ian,itemp1,nslt1,nslt2,0.0,spacing,corner,ecell,
     -  edge,xyzmin,xyzmax,vtot,nxyz,ixyz,itemp3,itemp2,indexs,ifirst,
     -  ilast,ntotcell,0,iout,maxrec)
      if (nconfig .eq. 1)
     -  print *,'Number of boxes in the X, Y, Z direction=',nxyz
      rmax=amax1(rsltmax,rcvmax)
      rmax2=rmax**2
      rcvmax2=rcvmax**2
      if (rsltmax .eq. 0.0) call zeroit(rnearsq,numsolv)
      call zeroit(cv,numsolv)
      call zeroiti(ianear,0,numsolv)
      incr=nslt
c     do is=1,1
      do is=1,numsolv
c       Set the grid limits for search
        call trnsfr(cslv,c(1,incr+iarepslv),3)
        ifar=0
        inside=0
        do k=1,3
          if (cslv(k) .lt. xyzmin(k)) then
            if (xyzmin(k)-cslv(k) .gt. rmax) then
              ifar=1
            else
              ixyz1(k)=0
              ixyz20=(cslv(k)+rmax-xyzmin(k))/ecell(k)
              ixyz2(k)=min0(nxyz(k)-1,ixyz20+1)
            end if
          else if (cslv(k) .ge. xyzmax(k)) then
            if (cslv(k)-xyzmax(k) .gt. rmax) then
              ifar=1
            else
              ixyz10=(cslv(k)-rmax-xyzmin(k))/ecell(k)
              ixyz1(k)=ixyz10-1
              ixyz2(k)=nxyz(k)-1
            end if
          else
            ixyz10=(cslv(k)-rmax-xyzmin(k))/ecell(k)
            ixyz1(k)=max0(0,ixyz10-1)
            ixyz20=(cslv(k)+rmax-xyzmin(k))/ecell(k)
            ixyz2(k)=min0(nxyz(k)-1,ixyz20+1)
            inside=inside+1
          end if
        end do
c       write (77,8791) is,ifar,cslv
c8791   format(i7,' IFAR=',i2,' CSLV=',3f10.5)
c       if (ifar .eq. 0) write (77,*) 'IXYZ1=',ixyz1
c       if (ifar .eq. 0) write (77,*) 'IXYZ2=',ixyz2
c       if (inside .eq. 3) then
c         Check the box enclosinf the solvent
c         ic=1+ix+iy*nxyz(1)+iz*nxyz(1)*nxyz(2)
c         if (ifirst(ic) .gt. 0) then
c           do ia=ifirst(ic),ilast(ic)
c             ix=(c(1,indexs(ia))-xyzmin(1))/ecell(1)
c             iy=(c(2,indexs(ia))-xyzmin(2))/ecell(2)
c             iz=(c(3,indexs(ia))-xyzmin(3))/ecell(3)
c             icc=1+ix+iy*nxyz(1)+iz*nxyz(1)*nxyz(2)
c             write (77,*) 'IXYZ(c),ic=',ix,iy,iz,icc
c           end do
c         else
c           write (77,*) 'EMPTY'
c         end if
c       end if
        rmin2=999999.0
        ncvsum=0
c       Use all cells for debug
c       do k=1,3
c         ixyz1(k)=0
c         ixyz2(k)=nxyz(k)-1
c       end do
c       ifar=0
        call zeroitd(cvsum,3)
        if (ifar .eq. 0) then
          do ix=ixyz1(1),ixyz2(1)
            do iy=ixyz1(2),ixyz2(2)
              do iz=ixyz1(3),ixyz2(3)
                ic=1+ix+iy*nxyz(1)+iz*nxyz(1)*nxyz(2)
                if (ifirst(ic) .gt. 0) then
                  do ia=ifirst(ic),ilast(ic)
                    rsvst2=dist2(c(1,indexs(ia)),c(1,incr+iarepslv))
                    if (rsvst2 .lt. rmin2) then
                      rmin2=rsvst2
                      ianear(is)=indexs(ia)
                    end if
                    if (rsvst2 .le. rmax2) then
                      if (cvmin .gt. 0.0) then
                        if (rsvst2 .le. rcvmax2) then
                          do k=1,3
                            d(k)=c(k,indexs(ia))-c(k,incr+iarepslv)
                          end do
c                         write (77,*) 'd=',d
                          call norm(d,1.0)
c                         write (77,*) 'dnorm=',d
                          do k=1,3
                            cvsum(k)=cvsum(k)+d(k)
                          end do
c                         write (77,*) 'cvsum=',cvsum
                          ncvsum=ncvsum+1
                        end if
                      end if
                    end if
                  end do
                end if
              end do
            end do
          end do
        end if
        incr=incr+naslv
        rnearsq(is)=rmin2
        if (ncvsum .gt. 0) cv(is)=1.d0-
     -    dsqrt(cvsum(1)**2+cvsum(2)**2+cvsum(3)**2)/dfloat(ncvsum)
        ncvsums(is)=ncvsum
        call trnsfrd(cvsums(1,is),cvsum,3)
        if (cv(is) .gt. 1.001) then
          print *,'PROGRAM ERROR: cv > 1.0:',cv(is),' ncvsum=',ncvsum
          stop
        end if
c       write (77,8945) is,sqrt(rnearsq(is)),cv(is),ncvsum
c8945   format(i8,' rnear=',f8.2,' cv=',f6.3,' ncvsm=',i7)
      end do
      return
      end
