      subroutine compact(co,cn,iatnum,ih,n,nnh,c0,rmin,rorgext,rorgcom,
     -  list,np,nlist,nocom,r,maxrec)
      dimension co(3,n),cn(3,n),iatnum(n),ih(n),c0(3),list(n),r(n)
c     Calculate the center and radius of the smallest sphere enclosing
c     all atoms in c
      character*2 iatnm2
      common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
     -  mmatno(64),iatnm2(99)
      dimension cmin(3),cmax(3),del(3),c1(3),c2(3),c3(3),c4(3),ca(3),
     -  rd(4),czero(3),cai(2),com(3),iopp(3)
      data iopp/0,0,0/,idiagplane /0/,nopp /0/
      call zeroit(czero,3)
c     print *,'COMPACT n,nnh=',n,nnh
      call extension(co,ih,nnh,1,n,cmin,cmax,c0,0,0,v)
      r0=0.0
      do k=1,3
c       print *,'c0,cmax,cmin=',c0(k),cmax(k),cmin(k)
        if (cmax(k) -c0(k) .gt. r0) r0=cmax(k)-c0(k)
      end do
      if (nocom .eq. 0) then
        do i=1,n
          r(i)=0.0
          if (iatnum(i) .gt. 0) r(i)=aw(iatnum(i))
        end do
        call cofms(co,com,n,r)
      end if
c     Calculate distances from c0, collect those within 60% of r0
      r02=r0**2
c     print *,'R02=',r02
      extfac=0.36
31    nlist=0
      rmaxext=0
      do ii=1,nnh
        i=ih(ii)
        ri2=dist2(co(1,i),c0)
c       print *,'i=',i,' ri2,r02=',ri2,r02
        if (ri2 .ge. r02*extfac) then
          nlist=nlist+1
          list(nlist)=i
          r(nlist)=ri2
          if (ri2 .gt. rmaxext) then
            rmaxext=ri2
            lmax=nlist
          end if
        end if
      end do
      if (nlist .lt. 2) then
        print *,'Too few atoms remained:',nlist
        call askyn('Do you want to keep all atoms instead',37,
     -    1,1,ikeep,0,0)
        if (ikeep .eq. 0) stop
        extfac=0.0
        go to 31
      end if
      nlstorg=nlist
      rorgext=sqrt(rmaxext)
      if (rorgext .lt. 0.01) then
        print *,'ERROR: The original radius of the solute is < 0.01'
        stop
      end if
c     print *,'nlist,rorgext=',nlist,rorgext
      if (nocom .eq. 0) then
        rmaxcom=0.0
        do ii=1,nnh
          i=ih(ii)
          ri2=dist2(co(1,i),com)
          if (ri2 .gt. rmaxcom) then
            rmaxcom=ri2
          end if
        end do
        rorgcom=sqrt(rmaxcom)
      end if
c     write (6,2344) (list(i),i=1,6),(r(i),i=1,6)
      nstep=1
      np=1
      call swapl(list,r,nlist,lmax,-1,maxrec)
      call lpshiftc(co,list,nlist,del,c0,np,maxrec)
c     write (6,2344) (list(i),i=1,6),(r(i),i=1,6)
      call findshift(co,r,list,c0,del,nlist+1,nlist,rlambda,lmax,nzero,
     -  maxrec)
c     print *,'nzero(1)=',nzero
c     print *,'lambda=',rlambda
      call transc(c0,del,rlambda)
      call newdist(r,co,list,c0,nlstorg,maxrec)
c     write (6,2344) (list(i),i=1,6),(r(i),i=1,6)
      np=2
      nlistdel=0
      call swapl(list,r,nlist,lmax,-1,maxrec)
30    nstep=nstep+1
      call lpshiftc(co,list,nlist,del,c0,np,maxrec)
c2344 format(' Nlist=',6i2,' r=',6f8.3)
      call findshift(co,r,list,c0,del,nlist+1,nlist-nlistdel,
     -  rlambda,lmax,nzero,maxrec)
c     print *,'nzero(2)=',nzero
c     print *,'lambda=',rlambda
      nlistdel=0
      do k=1,3
        c1(k)=co(k,list(nlist+1))-c0(k)
      end do
      rlammax=scprod(c1,del)
      if ((rlambda .ge. rlammax .or. rlambda .eq. 0.0) .and.
     -  nzero .le. 1) then
c       Make these two atoms the diameter
        write (6,1001) list(nlist+1),list(nlist+2)
        rlambda=rlammax
        go to 90
      end if
c     Check the triangle lmax, nlist+1, nlist+2 if it has an obtuse angle
33    call angles(r(lmax),r(nlist+1),r(nlist+2),cal,ca(1),ca(2))
      do ic=1,2
        if (ca(ic) .lt. 0.0) then
          call swapl(list,r,nlist+ic,lmax,0,maxrec)
          call swapl(list,r,nlist-nlistdel,lmax,0,maxrec)
          nlistdel=nlistdel+1
c         nlistdel>0 will make the lambda search bypass the last nlistdel ats
          call transc(c0,del,rlambda)
          call newdist(r,co,list,c0,nlstorg,maxrec)
          write (6,1002) list(nlist),list(nlist+1),list(nlist+2)
c         Go back to generate new triangle
          go to 30
        end if
      end do
      call transc(c0,del,rlambda)
      call newdist(r,co,list,c0,nlstorg,maxrec)
      np=3
c     write (6,6666) r(lmax),(r(nlist+ii),ii=1,np)
      if (nlistdel .eq. 0) then
        call swapl(list,r,nlist,lmax,-1,maxrec)
      else
        call swapl(list,r,nlist,nlist-nlistdel,0,maxrec)
        call swapl(list,r,nlist,lmax,-1,maxrec)
      end if
c     write (6,6666) r(lmax),(r(nlist+ii),ii=1,np)
c6666  format(' Distances of vertices from c0=',4f8.2)
40    nstep=nstep+1
c     Find shift preserving the distance from three points
      do k=1,3
        c1(k)=co(k,list(nlist+2))-co(k,list(nlist+1))
        c2(k)=co(k,list(nlist+3))-co(k,list(nlist+1))
      end do
      call vprd(c1,c2,del)
      do k=1,3
        c1(k)=co(k,list(nlist+1))-c0(k)
      end do
      facsgn=1.0
      if (scprod(c1,del) .lt. 0.0) facsgn=-1.0
      call norm(del,facsgn)
      call findshift(co,r,list,c0,del,nlist+1,nlist-nlistdel,
     -  rlambda,lmax,nzero,maxrec)
      nlistdel=0
c     Check for maximum
      rlammax=scprod(c1,del)
      idiagplane=0
      if (rlambda .ge. rlammax) then
c       First three atoms are in a diagonal plane
        write (6,1003) list(nlist+1),list(nlist+2),list(nlist+3)
        rlambda=rlammax
        idiagplane=1
        go to 90
      end if
c     Check for center of the sphere being outside the tetrahedron
      call transc(c0,del,rlambda)
      call newdist(r,co,list,c0,nlstorg,maxrec)
      rlambda=0.0
      ix=list(lmax)
      do ip=1,np
c       iopp(i) will be set to ip if atom (nlist+ip) and c0 are separated
c       by a tetrahedron face
        nopp=0
        call getindex(nlist,ip,iv,i1,i2)
        do k=1,3
          c1(k)=co(k,list(i1))-co(k,ix)
          c2(k)=co(k,list(i2))-co(k,ix)
          c3(k)=co(k,list(iv))-co(k,ix)
          c4(k)=c0(k)-co(k,ix)
        end do
        call vprd(c1,c2,ca)
        if (scprod(ca,c3)*scprod(ca,c4) .lt. 0.0) then
          nopp=nopp+1
          iopp(nopp)=ip
        end if
      end do
      if (nopp .eq. 1) then
c       First drop (temporarily) vertex nlist+iopp(1)
        ip=nlist+iopp(1)
c       swap vertex ip with nlist+1
        call swapl(list,r,nlist+1,ip,0,maxrec)
        write (6,1004) list(lmax),list(nlist+1),list(nlist+2),
     -    list(nlist+3)
c       rr=sqrt(r(nlist+1))
c       print *,'R=',rr
c       write (6,6666) r(lmax),(r(nlist+ii),ii=1,np)
        nlistdel=nlistdel+1
c       nlistdel>0 will make the lambda search bypass the last nlistdel ats
        nlist=nlist+1
        np=np-1
c       Go back to check the triangle
        go to 33
      else if (nopp .eq. 2) then
c       Two opposing vertices found
        do ic=1,2
          ip=iopp(ic)
          call getindex(nlist,ip,iv,i1,i2)
          call angles(dist2(co(1,i1),co(1,i2)),r(i1),r(i2),
     -                caix,cai(1),cai(2))
          if (cai(1)*cai(2) .gt. 0.0) then
c           Acute triangle face found, just drop ip
c           Swap vertex ip with lmax
            call swapl(list,r,nlist+ip,lmax,0,maxrec)
c           Move ip to the top of the list for temporary disregard
            call swapl(list,r,nlist,lmax,0,maxrec)
            nlistdel=nlistdel+1
            write (6,1005) list(nlist),
     -        list(nlist+1),list(nlist+2),list(nlist+3)
c           rr=sqrt(r(nlist+1))
c           print *,'R=',rr
c           write (6,6666) (r(nlist+ii),ii=1,np)
c           Go back to generate new tetrahedron
            go to 40
          end if
        end do
c       If reached here, neither faces are acute triangles - temporarily
c       drop both
        call swapl(list,r,nlist+iopp(1),nlist+1,0,maxrec)
        call swapl(list,r,nlist+iopp(2),nlist+2,0,maxrec)
        nlistdel=2
        nlist=nlist+2
        np=2
        go to 30
      end if
c     If reached here, center is inside the tetrahedron - done.
      np=4
      call swapl(list,r,nlist,lmax,-1,maxrec)
c     Calculation done, perform the shift
90    call transc(c0,del,rlambda)
      call newdist(r,co,list,c0,nlstorg,maxrec)
c     write (6,6666) (r(nlist+ii),ii=1,np)
      call shiftmol(co,n,c0,cn,-1.0)
c     Calculate the minimum radius (as the maximum of the surface atom radii)
      rminmax=0.0
      do ip=1,np
        ii=list(nlist+ip)
        rd(ip)=sqrt(cn(1,ii)**2+cn(2,ii)**2+cn(3,ii)**2)
        if (rminmax .lt. rd(ip)) rminmax=rd(ip)
      end do
      rmin=rminmax
      write (6,1000) (list(nlist+ip),ip=1,np)
      ndiff=0
      do ip=2,np
        if (abs(rd(ip)-rmin) .gt. 0.001) ndiff=ndiff+1
      end do
      if (ndiff-idiagplane .gt. 0) then
          print *,'PROGRAM ERROR: The radii of the points on',
     -    ' the sphere are different:'
         write (6,1010) (rd(ii),ii=1,np)
      end if
      call chksphr(cn,ih,n,nnh,rmin,czero)
      return
1000  format(' --- Sphere-spanning atoms:',4i6)
1001  format(' --- Diameter found spanned by atoms',2i6)
1002  format(' --- Center of circle outside the triangle',3i6)
1003  format(' --- Diagonal plane found spanned by atoms',3i6)
1004  format(' --- Sphere center is outside one tetrahedron face',4i6)
1005  format(' --- Sphere center is outside two tetrahedron faces',4i6)
1010  format('Radii=',4f8.4)
      end
