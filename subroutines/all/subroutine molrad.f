      subroutine molrad(c,index,nats,iout,maxats)
      dimension c(3,maxats),index(nats)
      parameter (MAXFRAMES=50000,MAXCOPY=600)
      common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
     -  increment,increment2,res(2,MAXFRAMES,MAXCOPY),
     -  x0(MAXCOPY),y0(MAXCOPY),nxselres,ixselres(MAXCOPY)
      real*8 rg,rh,rav(3),rgten(3,3),diag(3),offdiag(3),diagmin,diagmid,
     -  diagmax
      write (iout,*)
      call zeroitd(rav,3)
      call zeroitd(rgten,9)
      do ia=1,nats
        do k=1,3
          rav(k)=rav(k)+c(k,index(ia))
        end do
      end do
      do k=1,3
        rav(k)=rav(k)/nats
      end do
      rg=0.d0
      do ia=1,nats
        rg=rg+(c(1,index(ia))-rav(1))**2+(c(2,index(ia))-rav(2))**2+
     -    (c(3,index(ia))-rav(3))**2
      end do
      rg=rg/nats
      rg=sqrt(rg)
      rh=0.d0
      do ia=1,nats
        do ja=ia+1,nats
          rij=dist2(c(1,index(ia)),c(1,index(ja)))
          if (rij .lt. 1.e-5) then
            write (6,1000) ia,ja
          else
            rh=rh+1.0/sqrt(rij)
          end if
          do k=1,3
            do l=1,3
              rgten(k,l)=rgten(k,l)+(c(k,index(ia))-c(k,index(ja)))*
     -          (c(l,index(ia))-c(l,index(ja)))
            end do
          end do
        end do
      end do
      rh=2.0*rh/nats**2
c     Find the eigenvectors a and eigenvalues mu of rgten
      do k=1,3
        do l=1,3
          rgten(k,l)=rgten(k,l)/float(nats**2)
        end do
      end do
      call dtred2(rgten,3,3,diag,offdiag)
      call dtqli(diag,offdiag,3,3,rgten,ierr)
      if (ierr .gt. 0) then
        write (6,1004)
        write (iout,1004)
        return
      end if
      diagmin=dmin1(diag(1),diag(2),diag(3))
      diagmax=dmax1(diag(1),diag(2),diag(3))
      diagmid=0.d0
      if (diagmax .eq. diagmin) then
        diagmid=diagmin
      else
        do k=1,3
          if (diag(k) .ne. diagmin .and. diag(k) .ne. diagmax)
     -      diagmid=diag(k)
        end do
      end if
      rmaxmin=1.0
      if (diagmin .ne. 0.d0) rmaxmin=diagmax/diagmin
      asph=diagmax-0.5*(diagmin+diagmid)
      acyl=diagmid-diagmin
      rsa=1.5*(diag(1)**2+diag(2)**2+diag(3)**2)/
     -  (diag(1)+diag(2)+diag(3))**2-0.5
      if (nframe .gt. 0) then
        call trajlimtest(nframe,MAXFRAMES)
        write (iout,1001) nframe,rg,1/rh
        res(1,nframe,1)=rg
        res(2,nframe,1)=1.0/rh
        res(1,nframe,2)=rh
        res(2,nframe,2)=rmaxmin
      else
        write (iout,1002) rg,1/rh
      end if
      write (iout,1003) diag,rmaxmin,asph,acyl,rsa
      return
1000  format(' Atoms ',i6,' and ',i6,' are too close')
1001  format(' Frame',i6,' R(gyration)=',f9.3,' R(hydrodynamic)=',f9.3)
1002  format(' R(gyration)=',f9.3,' R(hydrodynamic)=',f9.3)
1003  format(' Moments of inertia=',3f10.3,' M(max)/M(min)=',f10.3,/,
     -  ' Asphericity=',f8.3,' Acylindricity=',f8.3,
     -  ' Relative shape anisotropy=',f8.3)
1004  format(' Calculation aborted due to diagonalization failure')
      end
