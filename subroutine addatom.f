      subroutine addatom(ityp,r4,r3,r2,r1,x,rij,aijk,tijkl,bend,ian,pi,
     -  iconv,ifail)
c     Add one atom with known distance, angle and torsion
      dimension r4(3),r3(3),r2(3),r1(3),x(3)
      dimension c12(3),c13(3),c14(3),c23(3),x0(3),e2(3),e3(3),bs(3)
      dimension p(3,4)
      character*13 lab
c     print *,'ADDATOM ityp=',ityp
      if (ian .gt. 0) then
        write (lab,1004) ian
        llab=13
      else
        lab=' X'
        llab=2
      end if
      if (iconv .eq. 1) then
        angconv=pi/180.0
        wconv=1.0
      else
        angconv=1.0
        wconv=180.0/pi
      end if
      if (ityp .eq. 1) then
c       Seqential definition
        write (6,1001) lab(1:llab),rij,wconv*aijk,wconv*tijkl
        call arrdiff(r1,r2,c12,3)
        call arrdiff(r3,r2,c23,3)
c       print *,'ADDATOM r1=',r1
c       print *,'ADDATOM r2=',r2
c       print *,'ADDATOM r3=',r3
        rnorm12=sqrt(scprod(c12,c12))
        rnorm23=sqrt(scprod(c23,c23))
        if (aijk .lt. 179.9) then
          if (abs(scprod(c12,c23)/(rnorm12*rnorm23)) .gt. 0.999) then
            print *,'Root atoms are colinear - new atom site is ',
     -        'undefined'
            ifail=ifail+1
          end if
          do k=1,3
            c12(k)=c12(k)/rnorm12
            c23(k)=c23(k)/rnorm23
            x0(k)=r1(k)-rij*c12(k)*cos(angconv*aijk)
          end do
          xproj=rij*sin(angconv*aijk)
          call vprd(c12,c23,e3)
          call vprd(e3,c12,e2)
          do k=1,3
            x(k)=x0(k)+xproj*(cos(angconv*tijkl)*e2(k)+
     -        sin(angconv*tijkl)*e3(k))
          end do
        else
c         Atoms r3, r2, r1, and x are colinear
          do k=1,3
            x(k)=r1(k)+rij*c12(k)/rnorm12
          end do
        end if
c       print *,'ADDATOM x=',x
        call trnsfr(p(1,1),r3,3)
        call trnsfr(p(1,2),r2,3)
        call trnsfr(p(1,3),r1,3)
        call trnsfr(p(1,4),x,3)
c       print *,'2,3,4 angle=',angleijk(p,4,2,3,4,6)*180.0/3.1415
c       print *,'1,2,3,4 dih angle=',
c    -    dihangl(p,1,2,3,4,0,100000)*180.0/3.1415
      else if (ityp .eq. 2) then
c       Bisector definition
c       print *,'BEND=',bend,' AIJK,TIJKL=',aijk,tijkl
        write (6,1002) lab(1:llab),rij
c       print *,'R1=',r1
c       print *,'R2=',r2
c       print *,'R3=',r3
c       print *,'RIJ=',rij
        call arrdiff(r1,r2,c12,3)
        call arrdiff(r1,r3,c13,3)
        rnorm12=sqrt(scprod(c12,c12))
        rnorm13=sqrt(scprod(c13,c13))
        do k=1,3
          c12(k)=c12(k)/rnorm12
          c13(k)=c13(k)/rnorm13
        end do
        colinear=scprod(c12,c13)
        if (colinear .lt. -0.999) then
          print *,'Atoms R2, R1, and R3 are colinear - can not use ',
     -      'bisector definition'
            ifail=ifail+1
        end if
        call arrsum(c12,c13,bs,3)
        rnormbs=sqrt(scprod(bs,bs))
        if (bend .eq. 0.0) then
          do k=1,3
            bs(k)=bs(k)*rij/rnormbs
            x(k)=r1(k)+bs(k)
          end do
        else
          call arrdiff(r2,r3,c23,3)
          call vprd(c23,bs,e3)
          rnorm3=sqrt(scprod(e3,e3))
          rnormbs=sqrt(scprod(bs,bs))
          bendr=bend*angconv
          do k=1,3
            x(k)=r1(k)+
     -        rij*(bs(k)*cos(bendr)/rnormbs+e3(k)*sin(bendr)/rnorm3)
          end do
        end if
c       print *,'X =',x
      else if (ityp .eq. 3) then
c       Trisector definition
        write (6,1003) lab(1:llab),rij
        call arrdiff(r1,r2,c12,3)
        call arrdiff(r1,r3,c13,3)
        call arrdiff(r1,r4,c14,3)
        rnorm12=sqrt(scprod(c12,c12))
        rnorm13=sqrt(scprod(c13,c13))
        rnorm14=sqrt(scprod(c14,c14))
        do k=1,3
          c12(k)=c12(k)/rnorm12
          c13(k)=c13(k)/rnorm13
          c14(k)=c14(k)/rnorm14
          bs(k)=c12(k)+c13(k)+c14(k)
        end do
        rnormbs=sqrt(scprod(bs,bs))
        if (rnormbs .lt. 0.01) then
          print *,'Atom R1 is in the plane of R2, R3, and R4 ',
     -      'can not use trisector definition'
            ifail=ifail+1
        end if
        do k=1,3
          bs(k)=bs(k)*rij/rnormbs
          x(k)=r1(k)+bs(k)
        end do
      end if
      return
1001  format(' Adding',a,': r(X-R1)=',f7.4,' a(R2-R1-X)=',f7.2,
     -  ' t(R3-R2-R1-X)=',f7.2)
1002  format(' Adding',a,': r(X-R1)=',f7.4,' along the negative ',
     -  'bisector of the R2-R1-R3 angle')
1003  format(' Adding',a,': r(X-R1)=',f7.4,' along the negative ',
     -  'trisector of the',/,' pyramid with base R2, R3, and R3')
1004  format(' atom #',i6)
      end
