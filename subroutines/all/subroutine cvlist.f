      subroutine cvlist(c,n,nslt,nsltref_f,nsltref_l,naslv,islvrep,
     -  icvtyp,rcut,rprox,cv,ixtemp,rtemp,dtemp,ifa,ila,it1,ctemp,
     -  isortslv,cvlim)
      dimension c(3,n),cv(n),ixtemp(n),rtemp(n),ifa(n),ila(n),it1(n),
     -  ctemp(3,n),dtemp(n),rprox(n)
      real*8 xnum,ynum,znum,den
      rcut2=rcut*rcut
c     Calculating the circular variance for solute atoms and solvent molecules
c     print *,'CVLIST n,nslt,naslv=',n,nslt,naslv
      call zeroit(rprox,nslt)
      do ia=1,nslt
        xnum=0.d0
        ynum=0.d0
        znum=0.d0
        den=0.d0
        do ja=nsltref_f,nsltref_l
          if (ia .ne. ja) then
            dx=c(1,ia)-c(1,ja)
            dy=c(2,ia)-c(2,ja)
            dz=c(3,ia)-c(3,ja)
            d2=dx*dx+dy*dy+dz*dz
            if (d2 .le. rcut2) then
              d=sqrt(d2)
              if (icvtyp .eq. 1) then
                xnum=xnum+dx/d
                ynum=ynum+dy/d
                znum=znum+dz/d
                den=den+1.d0
              else
                xnum=xnum+dx
                ynum=ynum+dy
                znum=znum+dz
                den=den+d
              end if
            end if
          end if
        end do
        if (den .gt. 0.d0) then
          cv(ia)=1.d0-dsqrt(xnum*xnum+ynum*ynum+znum*znum)/den
        else
          cv(ia)=0.0
        end if
      end do
      if (n .gt. nslt) then
c       CV for solvents
        numslv=(n-nslt)/naslv
        iaa=nslt+islvrep
        do ia=1,numslv
          xnum=0.d0
          ynum=0.d0
          znum=0.d0
          den=0.d0
          d2min=1000000.0
          do ja=nsltref_f,nsltref_l
            dx=c(1,iaa)-c(1,ja)
            dy=c(2,iaa)-c(2,ja)
            dz=c(3,iaa)-c(3,ja)
            d2=dx*dx+dy*dy+dz*dz
            if (d2 .le. rcut2) then
              d=sqrt(d2)
              if (icvtyp .eq. 1) then
                xnum=xnum+dx/d
                ynum=ynum+dy/d
                znum=znum+dz/d
                den=den+1.d0
              else
                xnum=xnum+dx
                ynum=ynum+dy
                znum=znum+dz
                den=den+d
              end if
            end if
            if (d2min .gt. d2) d2min=d2
          end do
          cvi=0.0
          if (den .gt. 0.d0)
     -      cvi=1.d0-dsqrt(xnum*xnum+ynum*ynum+znum*znum)/den
          d2min=sqrt(d2min)
          do i=nslt+(ia-1)*naslv+1,nslt+ia*naslv
            cv(i)=cvi
            rprox(i)=d2min
          end do
          iaa=iaa+naslv
          rtemp(ia)=cvi
          dtemp(ia)=sqrt(d2min)
          ixtemp(ia)=ia
        end do
        if (isortslv .eq. 1) then
          call trnsfr(ctemp(1,nslt+1),c(1,nslt+1),3*(n-nslt))
          call mrgsrt(6,ixtemp,rtemp,numslv,ifa,ila,it1,ctemp,n)
          limfound=0
          if (cvlim .ge. 1.0) limfound=1
          do ia=1,numslv
            ia0=nslt+(ia-1)*naslv
            ix0=nslt+(ixtemp(ia)-1)*naslv
            do i=ia0+1,ia0+naslv
              cv(i)=rtemp(ia)
              rprox(i)=dtemp(ia)
              ix=ix0+(i-ia0)
              call trnsfr(c(1,i),ctemp(1,ix),3)
            end do
            if (limfound .eq. 0) then
              if (rtemp(ia) .gt. cvlim) then
                limfound=1
                write (6,1000) cvlim,ia-1,cvlim,numslv-ia+1,numslv
              end if
            end if
          end do
        end if
      end if
      return
1000  format(' Number of solvents with CV < ',f6.3,'=',i4,/,
     -  ' Number of solvents with CV > ',f6.3,'=',i4,
     -  ' (total: ',i6,')')
      end
