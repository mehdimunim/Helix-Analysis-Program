      subroutine filterslv(c,ian,nslt,n,naslv,numsolv,molsltlim,nsegslt,
     -  index,nconfig,ifilttyp,intftyp,rsltmax,rcvmax,cvmin,ang12min1,
     -  ang12min2,r12max,cvrlim,rneigh,mol11,mol12,mol21,mol22,spacing,
     -  ixdrop,iarepslv,numsolvleft,nhbneig,ineig,ixres,nrecdel,ierr,
     -  iout,mxrec,maxng)
      dimension c(3,n),ian(n),index(n),molsltlim(3,nsegslt),
     -  ixdrop(mxrec),nhbneig(mxrec),ineig(maxng,mxrec),
     -  ixres(mxrec)
      parameter(MAXREC=200000,MAXPHI=400,IFILL10=MAXPHI**3-26*MAXREC)
      real*8 cvsums,cvsums2
      common /nnwork/ cvsums(3,MAXREC),cvsums2(3,MAXREC),
     -  ncvsums(MAXREC),ncvsums2(MAXREC),rnearsq(MAXREC),ianear(MAXREC),
     -  cv(MAXREC),rnearsq2(MAXREC),ianear2(MAXREC),cv2(MAXREC),
     -  indexs(MAXREC),ifirst(MAXREC),ilast(MAXREC),
     -  itemp1(MAXREC),itemp2(MAXREC),itemp3(MAXREC),fill(IFILL10)
      dimension r1(3),r2(3)
      real*8 cvsum(3)
c     print *,'FILTERSLV NUMSOLV=',numsolv,' NSLT,N=',nslt,n,
c    -  'ifilttyp=',ifilttyp
      nrecdel=0
      if (n .le. nslt) then
        write (iout,2001)
        ierr=1
        return
      end if
c     print *,'FILTERSLV n,nslt,naslv=',n,nslt,naslv
      if (ifilttyp .gt. 1) rsltmax=amax1(rsltmax,rcvmax,10.0)
      rsltmaxsq=rsltmax**2
      rcvmaxsq=rcvmax**2
      r12maxsq=r12max**2
      rneigh2=rneigh**2
      rminmin=999999.0
      cang12min1=cos(ang12min1*3.141592/180.0)
      cang12min2=cos(ang12min2*3.141592/180.0)
      call zeroiti(ixdrop,0,numsolv)
c     print *,'FILTERSLV numsolv,ifilttyp=',numsolv,ifilttyp
c     print *,'FILTERSLV rsltmax,rcvmax=',rsltmax,rcvmax
      if (ifilttyp .le. 3) then
        call calc_cv_rmin(c,ian,nslt,1,nslt,n,numsolv,naslv,iarepslv,
     -    spacing,rsltmax,rcvmax,cvmin,rnearsq,ianear,cv,cvsums,ncvsums,
     -    nconfig,indexs,ifirst,ilast,itemp1,itemp2,itemp3,iout,mxrec)
        do is=1,numsolv
          if (ifilttyp .gt. 1 .and. cv(is) .le. cvmin) ixdrop(is)=1
          if (rsltmax .gt. 0.0 .and. rnearsq(is) .gt. rsltmaxsq)
     -        ixdrop(is)=1
          if (rnearsq(is) .lt. rminmin) rminmin=rnearsq(is)
c         write (78,8921) is,cv(is),rnearsq(is),ixdrop(is),rminmin
c8921     format(i7,' cv=',f8.4,' rnsq=',f10.3,' ixdrop=',i2,
c    -      ' rminmin=',f8.2)
        end do
      else if (ifilttyp .eq. 4) then
        do is=1,numsolv
          ixdrop(is)=1
        end do
        nkeep=0
        incr=nslt
        idebug=0
        mol2fst=mol21
        do mol1=mol11,mol12
          call calc_cv_rmin(c,ian,nslt,molsltlim(1,mol1),
     -      molsltlim(2,mol1),n,numsolv,naslv,iarepslv,spacing,rsltmax,
     -      rcvmax,cvmin,rnearsq,ianear,cv,cvsums,ncvsums,nconfig,
     -      indexs,ifirst,ilast,itemp1,itemp2,itemp3,iout,mxrec)
          if (intftyp .eq. 3) mol2fst=mol1+1
          do mol2=mol2fst,mol22
            if (mol1 .ne. mol2) then
              call calc_cv_rmin(c,ian,nslt,molsltlim(1,mol2),
     -          molsltlim(2,mol2),n,numsolv,naslv,iarepslv,spacing,
     -          rsltmax,rcvmax,cvmin,rnearsq2,ianear2,cv2,cvsums2,
     -          ncvsums2,nconfig,indexs,ifirst,ilast,itemp1,itemp2,
     -          itemp3,iout,mxrec)
              nkeepij=0
              do is=1,numsolv
                ikeep=0
                if (idebug .gt. 0) write (76,9671) is,rnearsq(is),
     -            rnearsq2(is),ianear(is),ianear2(is)
                if (ianear(is)*ianear2(is) .gt. 0) then
                  d12sq=dist2(c(1,ianear(is)),c(1,ianear2(is)))
                  if (idebug .gt. 0) write (76,*) is,' D=',d12sq
                  if (d12sq .le. r12maxsq) then
c                   Solute anchor atoms are closer than r12max (15)
                    rlpsum=0.0
                    do k=1,3
                      r1(k)=c(k,ianear(is))-c(k,incr+iarepslv)
                      r2(k)=c(k,ianear2(is))-c(k,incr+iarepslv)
                      rlpsum=rlpsum+r1(k)*r2(k)
                    end do
                    ac=rlpsum/sqrt(abs(rnearsq(is)*rnearsq2(is)))
                    r1sq=r1(1)**2+r1(2)**2+r1(3)**2
                    r2sq=r2(1)**2+r2(2)**2+r2(3)**2
                    if (idebug .gt. 0) write (76,*) is,' AC=',ac,
     -                ' R1SQ,R2SQ=',r1sq,r2sq
                    cang12min=cang12min1
                    if (d12sq .gt. 36.0) cang12min=cang12min2
                    if (ac .lt. cang12min) then
c                     Solute-protein line angle > ang12min
                      if (idebug .gt. 0)
     -                  write (76,9672) is,cv(is),cv2(is),amax1(cv(is),
     -                    cv2(is))/amin1(cv(is),cv2(is))
                      if (amax1(cv(is),cv2(is))/amin1(cv(is),cv2(is))
     -                    .lt. cvrlim) then
c                       cv ratio is < cvrlim (3.5)
c                       Calculate combined CV
                        do k=1,3
                          cvsum(k)=cvsums(k,is)+cvsums2(k,is)
                        end do
                        cv12=1.d0-
     -                    dsqrt(cvsum(1)**2+cvsum(2)**2+cvsum(3)**2)/
     -                    dfloat(ncvsums(is)+ncvsums2(is))
                        if (idebug .gt. 0) write (76,*) is,' CV12=',cv12
                        if (cv12 .gt. cvmin) then
c                         CV wrt both proteins is > cv12lim (0.6)
                          ikeep=1
                        end if
                      end if
                    end if
                  end if
                end if
                if (idebug .gt. 0) write (76,*) is,' IXDROP=',ixdrop(is)
                if (rnearsq(is) .lt. rminmin) rminmin=rnearsq(is)
                if (rnearsq2(is) .lt. rminmin) rminmin=rnearsq2(is)
                incr=incr+naslv
                if (ikeep .eq. 1) nkeepij=nkeepij+1
                if (ikeep .eq. 1 .and. ixdrop(is) .eq. 1) then
                  ixdrop(is)=0
                  nkeep=nkeep+1
                  itemp1(nkeep)=is
                end if
              end do
              if (intftyp .gt. 1)
     -          write (6,2005) nconfig,mol1,mol2,nkeepij
            end if
          end do
        end do
c       Check for empty neighborhood
        nempty=0
        do is=1,numsolv
          if (ixdrop(is) .eq. 0) then
            if (amax1(rnearsq(is),rnearsq2(is)) .gt. rneigh2) then
              nn=0
              do iss=1,nkeep
                is2=itemp1(iss)
                if (is .ne. is2 .and. ixdrop(iss) .eq. 0) then
                  if (dist2(c(1,nslt+(is-1)*naslv+iarepslv),
     -              c(1,nslt+(is2-1)*naslv+iarepslv)) .le. rneigh2)
     -                nn=nn+1
                end if
              end do
              if (nn .eq. 0) then
                ixdrop(is)=1
                nempty=nempty+1
              end if
            end if
          end if
        end do
c       do is=1,numsolv
c         Calculate combined CV
c         do k=1,3
c           cvsum(k)=cvsums(k,is)+cvsums2(k,is)
c         end do
c         cv12s=1.d0-dsqrt(cvsum(1)**2+cvsum(2)**2+cvsum(3)**2)/
c    -      dfloat(ncvsums(is)+ncvsums2(is))
c         write (6,9681) is,ixdrop(is),cv12s,
c    -      (c(k,nslt+(is-1)*naslv+iarepslv),k=1,3)
c9681     format(i4,' ixdrop=',i2,' cv=',f8.5,' c=',3f10.5)
c       end do
      else if (ifilttyp .eq. 5) then
c       Hydrogen-bond bridging water search
        ia0=nslt
        do is=1,numsolv
          nhb=0
          ixdrop(is)=1
          ixres1=0
          do iaa=1,naslv
            ia=ia0+iaa
            do ihb=1,nhbneig(ia)
              ihbslt=ineig(maxng-ihb+1,ia)
              if (ihbslt .lt. nslt) then
                if (nhb .eq. 0) then
                  nhb=nhb+1
                  ixres1=ixres(ihbslt)
                else if (ixres1 .ne. ixres(ihbslt)) then
c                 2nd Hbond was found to a different residue
                  ixdrop(is)=0
                  go to  100
                end if
              end if
            end do
          end do
100       ia0=ia0+naslv
        end do
      end if
      if (nempty .gt. 0 .and. ifilttyp .eq. 4)
     -   write (iout,2002) nconfig,nempty
      incr=nslt
c     print *,'INDEX_UPDATE incr=',incr,' ix1=',index(incr+1)
      ndrop=0
      do is=1,numsolv
        if (ixdrop(is) .eq. 1) then
          ndrop=ndrop+1
        else
          do j=1,naslv
            index(incr+j-naslv*ndrop)=index(incr+j)
          end do
        end if
        incr=incr+naslv
      end do
      nrecdel=nrecdel+ndrop*naslv
      numsolvleft=numsolv-ndrop
      if (ifilttyp .lt. 5) write (iout,2003) nconfig,sqrt(rminmin)
      if (ndrop .eq. 0) then
        write (iout,*) 'NOTE: no solvent was dropped'
      else
        if (ndrop .eq. numsolv) then
          write (iout,*) 'All solvents were dropped'
        else
          write (iout,2004) nconfig,ndrop,numsolvleft
        end if
      end if
      return
2001  format(' ERROR: there are no solvents to filter')
2002  format(' Nconfig=',i6,':',i5,' solvents had empty neighbor ',
     -  'sphere')
2003  format(' Nconfig=',i6,' Closest solute-solvent distance=',f6.2,
     -  ' A')
2004  format(' Nconfig=',i6,' Number of solvents dropped=',i7,
     -  ' kept=',i6)
2005  format(' Nconfig=',i6,' Number of interface solvents between ',
     -  'molecules',i4,' and',i4,'=',i3)
9671  format(i5,' RNSQ,2=',2f10.3,' INEAR,2=',2i8)
9672  format(i5,' cv,cv2=',2f6.4,' R=',f8.4)
      end
