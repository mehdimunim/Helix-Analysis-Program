      subroutine pseudorot(a,n,ix,nmem,nneig,ineig,psr5,psr,npsr,ansrun,
     -  incgen,z,r,sinpsrs,cospsrs,q,qs,q2s,zav,zsq,zavs,zsqs,
     -  iprint,iout,nconf,radtodeg,maxng,maxrec)
      real*8 sinpsrs,cospsrs,qs,q2s,zavs,zsqs
      dimension a(3,n),ix(nmem),nneig(n),ineig(maxng,n),psr5(5),
     -  psr(nmem),z(nmem),r(3,nmem),sinpsrs(nmem),cospsrs(nmem),
     -  q(nmem),qs(nmem),q2s(nmem)
      character*1 ansrun
      dimension r0(3),r1(3),r2(3),rnorm(3),theta(5)
c     Calculate the pseudorotation angle defined by atoms i1-inmem,
      if (nmem .eq. 0) return
      npsr=0
      pi=atan(1.0)*4.0
      denomfac=2.0*(sin(36.0*pi/180.0)+sin(72.0*pi/180.0))
c     Check if atoms ix(1-nmem) form a loop
      do i=1,nmem
        if (ix(i) .lt. 1 .or. ix(i) .gt. n) then
          print *,'ERROR: invalid atomindex:',ix(i)
          return
        end if
        j=mod(i,nmem)+1
        do in=1,nneig(ix(i))
          if (ineig(in,ix(i)) .eq. ix(j)) go to 100
        end do
        write (iout,1000) ix(i),ix(j),sqrt(dist2(a(1,ix(i)),a(1,ix(j))))
100     continue
      end do
      npsr=(nmem-1)/2
      if (nconf .eq. 1) then
        zsqs=0.d0
        zavs=0.d0
        do i=1,npsr
          sinpsrs(i)=0.d0
          cospsrs(i)=0.d0
          qs(i)=0.d0
          q2s(i)=0.d0
        end do
      end if
      if (nmem .eq. 5) then
c       As defined by Altona & Sundaralingam, JACS (1972) 94, 8205.
c       Average over all cyclic permutations
        do it=1,5
          theta(it)=dihangl(a,ix(mod(it-1,5)+1),ix(mod(it,5)+1),
     -      ix(mod(it+1,5)+1),ix(mod(it+2,5)+1),0,maxrec)
        end do
c        write (iout,8711) (theta(it)*radtodeg,it=1,5)
c8711    format(' theta 1-5=',5f8.3)
        do id=0,4
          psr0=radtodeg*atan(((theta(mod(2+id,5)+1)+
     -         theta(mod(4+id,5)+1))-(theta(mod(1+id,5)+1)+
     -         theta(mod(3+id,5)+1)))/(theta(mod(id,5)+1)*denomfac))
          if (theta(mod(id,5)+1) .lt. 0) psr0=psr0+180.0
          psr5(id+1)=psr0
        end do
      end if
c     Cremer & Pople, JACS 97, 1354 (1975)
      call zeroit(r0,3)
      call zeroit(r1,3)
      call zeroit(r2,3)
      do i=1,nmem
        do k=1,3
          r0(k)=r0(k)+a(k,ix(i))
        end do
      end do
      do i=1,nmem
        do k=1,3
          r(k,i)=a(k,ix(i))-r0(k)/nmem
          r1(k)=r1(k)+r(k,i)*sin(2.0*pi*(i-1)/float(nmem))
          r2(k)=r2(k)+r(k,i)*cos(2.0*pi*(i-1)/float(nmem))
        end do
      end do
      call vprd(r1,r2,rnorm)
      rnmag=sqrt(rnorm(1)**2+rnorm(2)**2+rnorm(3)**2)
      do k=1,3
        rnorm(k)=rnorm(k)/rnmag
      end do
c      write (iout,2000) rnorm
c2000  format(' rnorm=',3f10.6)
      zsq=0.0
      do i=1,nmem
        z(i)=scprod(rnorm,r(1,i))
        zsq=zsq+z(i)**2
      end do
      zav=sqrt(zsq)
      if (nconf .gt. 0) then
        zavs=zavs+zav
        zsqs=zsqs+zsq
      end if
      sinsum=0.0
      cossum=0.0
      do m=2,npsr
        do i=1,nmem
          sinsum=sinsum-z(i)*sin(2*pi*m*(i-1-incgen)/float(nmem))
          cossum=cossum+z(i)*cos(2*pi*m*(i-1-incgen)/float(nmem))
        end do
        q(m)=sqrt((2.0/float(nmem))*(sinsum**2+cossum**2))
        tn=sinsum/cossum
        sqdenom=sqrt(1.0+tn*tn)
        sn=tn/sqdenom
        cs=1.0/sqdenom
        if (sinsum .lt. 0 .and. sn .gt. 0.0) sn=-sn
        if (cossum .lt. 0 .and. cs .gt. 0.0) cs=-cs
        if (nconf .gt. 0) then
          sinpsrs(m)=sinpsrs(m)+sn
          cospsrs(m)=cospsrs(m)+cs
          qs(m)=qs(m)+q(m)
          q2s(m)=q2s(m)+q(m)**2
        end if
        psr(m)=atansc(sn,cs,1,radtodeg)
      end do
      if (npsr .gt. 0)  then
        if (nmem .eq. 5) then
          psrndb=psr5(1)
c         write (iout,1004) psr5
          do id=1,4
            psrcorr=psr5(id+1)-id*144.0
            if (abs(psrcorr+360.0-psr5(id+1)) .lt.
     -          abs(psrcorr-psr5(id+1))) psrcorr=psrcorr+360.0
            if (abs(psrcorr-360.0-psr5(id+1)) .lt.
     -          abs(psrcorr-psr5(id+1))) psrcorr=psrcorr-360.0
            psrav=psrav+psrcorr
            psr5(id+1)=psrcorr
          end do
          psrav=psrav/5.0
        end if
        if (iprint .eq. 1) then
          if (ansrun .eq. '5') write (iout,1003) psrndb
          if (nmem .eq. 1) write (iout,1005) psr5,psrav
          write (iout,1007) zav
          write (iout,1006) 'angles',(psr(i),i=2,npsr)
          write (iout,1006) 'amplitudes',(q(i),i=2,npsr)
        end if
      end if
      return
1000  format(' WARNING: atoms ',i6,' and ',i6,' are not bonded, dij=',
     -  f8.3,' A',/,10x,'- ring is not complete')
1003  format(' Torsion-based pseudorotation angle, NDB convention=',
     -  f8.3,' deg')
c1004  format(' Torsion-based pseudorotation angles with all origins:',
c     -  /,5x,5f8.3,' deg')
1005  format(' Phase corrected torsion-based pseudorotation angles ',
     -  'with all origins:',/,5x,5f8.3,' deg   Mean=',f8.3,' deg')
1006  format(' General puckering ',a,' (Cremer & Pople)=',/,(10f8.2))
1007  format(' Root mean square distance from the mean plane=',f8.2,
     -  ' A')
      end
