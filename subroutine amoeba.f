      subroutine amoeba(p,y,mp,np,ndim,ftol,iter,c,cnew,ih,n,nnh,edge,
     -  ioppbc,cell,ncell,ixyzhex,rot,mintyp)
C*****Numerical Recipes
      parameter (NMAX=20,ALPHA=1.0,BETA=0.5,GAMMA=2.0,ITMAX=200)
      dimension p(MP,NP),y(MP),pr(NMAX),prr(NMAX),pbar(NMAX)
      dimension c(3,n),cnew(3,n),ih(n),edge(3),cell(3,27),
     -  rot(3,3),ixyzhex(3)
c     print *,'--- Simplex optimization started'
      mpts=ndim+1
      iiter=0
110   iter=0
1     ilo=1
      if (y(1) .gt. y(2)) then
        ihi=1
        inhi=2
      else
        ihi=2
        inhi=1
      end if
      do i=1,mpts
        if (y(i) .lt. y(ilo)) ilo=i
        if (y(i) .gt. y(ihi)) then
          inhi=ihi
          ihi=i
        else if (y(i) .gt. y(inhi)) then
          if (i .ne. ihi) inhi=i
        end if
      end do
c     write (6,1711) y,((p(i,k),k=1,3),i=1,4)
c1711  format(' Y=',4f10.5,/,(' P=',3f10.6))
      rtol=2.0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))
      if (rtol .lt. ftol) then
        call euler(rot,p(ilo,1),p(ilo,2),p(ilo,3))
        return
      end if
      if (iter .eq. itmax) then
        iiter=iiter+1
        if (iiter .gt. 2) then
          print *,'Too many resets'
          call euler(rot,p(ilo,1),p(ilo,2),p(ilo,3))
          return
        end if
        print *,'Reset ',iiter
        ybest=y(ilo)
        do i=1,3
          pr(i)=p(ilo,i)+i*0.05
          p(ihi,i)=pr(i)
        end do
        y(ihi)=touch(c,ih,n,nnh,edge,ioppbc,cell,ncell,ixyzhex,PR,
     -    rot,mintyp,0,cnew)
c       print *,'Reset ',iiter,' ilo=',ilo,' ihi=',ihi
        go to 110
      end if
      iter=iter+1
      ylo=-y(ilo)
      if (mod(iter,10) .eq. 0) write (6,1010) iter,ylo,ftol
1010  format(' --- Iteration',i4,' Objective function=',e15.8,' A',
     -  ' (tolerance=',f9.7,')')
      do j=1,ndim
        pbar(j)=0.0
      end do
      do i=1,mpts
        if (i .ne. ihi) then
          do j=1,ndim
            pbar(j)=pbar(j)+p(i,j)
          end do
        end if
      end do
      do j=1,ndim
        pbar(j)=pbar(j)/ndim
        pr(j)=(1.0+alpha)*pbar(j)-alpha*p(ihi,j)
      end do
c     ypr=funk(pr)
      YPR=touch(c,ih,n,nnh,edge,ioppbc,cell,ncell,ixyzhex,PR,rot,
     -  mintyp,0,cnew)
      if (ypr .le. y(ilo)) then
        do j=1,ndim
          prr(j)=gamma*pr(j)+(1.0-gamma)*pbar(j)
        end do
c       yprr=funk(prr)
        yprr=touch(c,ih,n,nnh,edge,ioppbc,cell,ncell,ixyzhex,PRR,rot,
     -  mintyp,0,cnew)
        if (yprr .lt. y(ilo)) then
          do j=1,ndim
            p(ihi,j)=prr(j)
          end do
          y(ihi)=yprr
        else
          do j=1,ndim
            p(ihi,j)=pr(j)
          end do
          y(ihi)=ypr
        end if
      else if (ypr .ge. y(inhi)) then
        if (ypr .lt. y(ihi)) then
          do j=1,ndim
            p(ihi,j)=pr(j)
          end do
          y(ihi)=ypr
        end if
        do j=1,ndim
          prr(j)=beta*p(ihi,j)+(1.-beta)*pbar(j)
        end do
c       yprr=funk(prr)
        yprr=touch(c,ih,n,nnh,edge,ioppbc,cell,ncell,ixyzhex,PRR,
     -    rot,mintyp,0,cnew)
        if (yprr .lt. y(ihi)) then
          do j=1,ndim
            p(ihi,j)=prr(j)
          end do
          y(ihi)=yprr
        else
          ywrst=y(ihi)
          ncont=0
240       continue
          ncont=ncont+1
c         write (6,1711) y,((p(i,k),k=1,3),i=1,4)
          yw=-10000.0
          do i=1,mpts
            if (i .ne. ilo) then
              do j=1,ndim
                pr(j)=0.5*(p(i,j)+p(ilo,j))
                p(i,j)=pr(j)
              end do
c             y(i)=funk(pr)
              y(i)=touch(c,ih,n,nnh,edge,ioppbc,cell,ncell,ixyzhex,PR,
     -          rot,mintyp,0,cnew)
              if (y(i) .gt. yw) yw=y(i)
            endiF
          end do
c         Continue contracting if the newest worst point is higher than
c         the old worst point
c         print *,'yw=',yw,' ywrst=',ywrst
          if (yw-ywrst .gt. ftol .and. ncont .lt. 5) go to 240
        end if
      else
        do j=1,ndim
          p(ihi,j)=pr(j)
        end do
        y(ihi)=ypr
      end if
      go to 1
      end
