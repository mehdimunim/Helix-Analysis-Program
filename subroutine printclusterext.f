      subroutine printclusterext(c,cent,ncl,ifclst,ilclst,ixclst,
     -  na,iclstyp,outfile,namleno,iout)
      dimension c(3,na),cent(3,ncl),ifclst(ncl),ilclst(ncl),ixclst(na)
      character*(*) outfile
      dimension xyzmin(3),xyzmax(3),xyz(3)
      character*1 xyzl
      character*200 centfile
      common /axislab/ xyzl(3)
      real*8 dxyz(3),r,theta,phi
c     print *,'PRINTCLUSTEREXT ncl=',ncl,' na=',na,
c    -  ' iclstyp,iout=',iclstyp,iout
      do icl=1,ncl
        do k=1,3
          xyzmin(k)=99999.0
          xyzmax(k)=-99999.0
        end do
        rmin=99999.0
        rmax=0.0
        thetamin=10.0
        thetamax=-10.0
        phimin=10.0
        phimax=-10.0
        do ia=ifclst(icl),ilclst(icl)
          do k=1,3
            xyz(k)=c(k,ixclst(ia))
            dxyz(k)=xyz(k)
            if (xyzmin(k) .gt. xyz(k)) xyzmin(k)=xyz(k)
            if (xyzmax(k) .lt. xyz(k)) xyzmax(k)=xyz(k)
          end do
c         Convert xyz to polar coordinates
          call polar(dxyz(1),dxyz(2),dxyz(3),r,phi,theta)
          if (r .lt. rmin) rmin=r
          if (r .gt. rmax) rmax=r
          if (theta .lt. thetamin) thetamin=theta
          if (theta .gt. thetamax) thetamax=theta
          if (phi .lt. phimin) phimin=phi
          if (phi .gt. phimax) phimax=phi
        end do
        write (iout,1000) icl
        if (iclstyp .ne. 2 .and. iclstyp .ne. 4) then
c         Calculate cluster COM
          call zeroit(cent(1,icl),3)
          do ia=ifclst(icl),ilclst(icl)
            do k=1,3
              cent(k,icl)=cent(k,icl)+c(k,ixclst(ia))
            end do
          end do
          do k=1,3
            cent(k,icl)=cent(k,icl)/float(ilclst(icl)-ifclst(icl)+1)
          end do
        end if
        write (iout,1002) (cent(k,icl),k=1,3)
        write (iout,1001) ('    '//xyzl(k),xyzmin(k),xyzmax(k),'A',
     -    k=1,3)
        write (iout,1001) '    R',rmin,rmax,'A'
        write (iout,1001) 'Theta',thetamin,thetamax,'rad'
        write (iout,1001) '  Phi',phimin,phimax,'rad'
      end do
      call askyn('Do you want a cluster-center file',33,
     -  1,-1,icentfile,000,0)
      if (icentfile .eq. 1) then
        centfile(1:namleno)=outfile(1:namleno)
        centfile(namleno+1:namleno+4)='.cnt'
        print *,'Center coordinates will be written to file ',
     -    centfile(1:namleno+4)
        call openfile(51,0,'cluster center',14,'new',centfile,namleno+4,
     -    notfound,0,1,1,0,0)
        write (51,1003) ncl,((cent(k,icl),k=1,3),icl=1,ncl)
      end if
      return
1000  format(/,' Cluster ',i5,':')
1001  format(1x,a,'-range: [',f6.1,' - ',f6.1,'] ',a)
1002  format(' Cluster center coordinates: ',3f10.5)
1003  format(i3,/,(3f12.5))
      end
