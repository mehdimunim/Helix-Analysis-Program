      subroutine clstrs_density(nfrst,n,iclst,ifirst,ilast,nofcls,
     -  nofcls_t,nng,ing,nng2r,idenclstyp,nnmin,it1,it2,it3,it4,it5,
     -  it6,mx2d,iout)
c*****Find all clusters in a network and sort atoms in a cluster by groups
      dimension iclst(mx2d),ifirst(mx2d),ilast(mx2d),nng(mx2d),
     -  ing(mx2d,mx2d),nng2r(mx2d),it1(mx2d),it2(mx2d),it3(mx2d),
     -  it4(mx2d),it5(mx2d),it6(mx2d)
c     Input parameters:
c     n0,n: Use atomnumbers (vertices) from n0 to n (inclusive)
c     idenclst=1: bond if share nnmin neighbors
c     idenclst=2: bond is within cutoff and share at least nnmin neighbors
c     idenclst=3: bond if within cutoff and each has at least nnmin neighbors
c     Workspace arrays: it1,it2,it3,it4,it5,it6
c     Output parameters
c     nofcls: Number of clusters requested
c     iclst,ifirst,ilast: The elements of the ig-th group are atoms
c     (iclst(ia),ia=ifirst(ig),ilast(ig))
c     Initialization
      print *
c     print *,'CLSTRS_DENSITY nnmin=',nnmin
      iverb=1
c     Find neighbor number range
      nngmin=mx2d
      nngmax=0
      nz=0
      do i=nfrst,n
        if (nng(i) .lt. nngmin) nngmin=nng(i)
        if (nng(i) .gt. nngmax) nngmax=nng(i)
        if (nng(i) .lt. nnmin) nz=nz+1
      end do
      write (iout,2001) nngmin,nngmax
      write (   6,2001) nngmin,nngmax
      if (nz .gt. 0) write (iout,2000) nnmin,nz
      if (idenclstyp .eq. 3) then
        do i=nfrst,n
          do jj=1,nng(i)
            j=iabs(ing(jj,i))
            if (i .lt. j) then
              if (nng(i) .lt. nnmin .or. nng(j) .lt. nnmin) then
c               Remove i-j bond
                do ii=1,nng(j)
                  if (ing(ii,j) .eq. i) ing(ii,j)=-ing(ii,j)
                end do
                ing(jj,i)=-ing(jj,i)
              end if
            end if
          end do
        end do
c        do i=nfrst,n
c          write (77,7943) i,(ing(j,i),j=1,nng(i))
c7943      format(i5,' NN:',50i5)
c        end do
        call clean_ng(nfrst,n,nng,ing,mx2d)
c        print *,'CLEAN done'
c        do i=nfrst,n
c          write (78,7943) i,(ing(j,i),j=1,nng(i))
c        end do
        call checknnlist(nfrst,n,ing,nng,nerr,mx2d)
        if (nerr .gt. 0) stop
        iverb=0
        call clstrs(ing,nng,it1,nfrst,n,iclst,ifirst,ilast,0,nofcls,
     -    it2,0,it2,0,0006,inperr,iverb,n,n,mx2d)
c       print *,'CLSTRS done'
      else if (idenclstyp .eq. 2) then
c       Remove bonds if the # of common neighbors is < nnmin
        do i=nfrst,n
          do j=1,nng(i)
            nmatch=0
            do ii=1,nng(i)
              do jj=1,nng(j)
                iii=iabs(ing(jj,j))
                if (iii .gt. i) then
                  if (iabs(ing(ii,i)) .eq. iii) nmatch=nmatch+1
                end if
              end do
            end do
            if (nmatch .lt. nnmin) then
c             Mark bond for removal by setting it to -ing(ii,i)
              do ii=1,nng(i)
                if (ing(ii,i) .eq. j) ing(ii,i)=-ing(ii,i)
              end do
              do jj=1,nng(j)
                if (ing(jj,j) .eq. i) ing(jj,j)=-ing(jj,j)
              end do
            end if
          end do
        end do
c       Remove the negative and zero ing entries
        call clean_ng(nfrst,n,nng,ing,mx2d)
        call clstrs(ing,nng,it1,nfrst,n,iclst,ifirst,ilast,0,nofcls,
     -    it2,0,it2,0,iout,inperr,0,n,n,mx2d)
      else if (idenclstyp .eq. 1) then
c       Keep bond when # of common neighbors is >= nnmin; use nng2r list too
        call zeroiti(it3,nfrst-1,n)
c       Move the R-2R neighbors together with the R neighbors
        do i=nfrst,n
          do j=1,nng2r(i)
            ing(nng(i)+j,i)=ing(mx2d-j+1,i)
          end do
        end do
        do i=nfrst,n
          do j=1,nng(i)+nng2r(i)
            jn=ing(j,i)
c           Check for # of common neighbors
            nmatch=0
            do ii=1,nng(i)
              do jj=1,nng(jn)
                if (ing(ii,i) .eq. ing(jj,jn)) nmatch=nmatch+1
              end do
            end do
            if (nmatch .ge. nnmin) then
c             Create bond
              ing(mx2d-it3(i),i)=jn
              ing(mx2d-it3(jn),jn)=i
              it3(i)=it3(i)+1
              it3(jn)=it3(jn)+1
            end if
          end do
        end do
c       Move the it3 bonds to the start in ing
        do i=nfrst,n
          do j=1,it3(i)
            ing(j,i)=ing(mx2d-j+1,i)
          end do
        end do
        call trnsfi(nng(nfrst),it3(nfrst),n-nfrst+1)
        call clstrs(ing,nng,it1,nfrst,n,iclst,ifirst,ilast,0,nofcls,
     -    it2,0,it2,0,iout,inperr,0,n,n,mx2d)
      end if
c     Sort clusters by the # of members
      do i=1,nofcls
        it6(i)=float(ifirst(i)-ilast(i))
      end do
      call indexit(it1,1,nofcls,0)
      call mrgsrti(iout,it1,it6,nofcls,it2,it3,it4,it5,nofcls)
c     print *,'MRGSRT done'
c      write (6,8943) (it1(i),i=1,nofcls)
c8943  format(' IT1:',/,(20i4))
c      write (6,8944) (-t2(i),i=1,nofcls)
c8944  format(' T2:',/,(20f4.0))
      inc=0
      do i=1,nofcls
        nmem=ilast(it1(i))-ifirst(it1(i))+1
        do j=1,nmem
          it4(inc+j)=iclst(ifirst(it1(i))-1+j)
        end do
        it2(i)=inc+1
        it3(i)=inc+nmem
        inc=inc+nmem
      end do
      call trnsfi(ifirst,it2,nofcls)
      call trnsfi(ilast,it3,nofcls)
      call trnsfi(iclst,it4,n-nfrst+1)
c      do i=1,nofcls
c        write (79,7891) i,ifirst(i),ilast(i),
c     -    (iclst(j),j=ifirst(i),ilast(i))
c7891    format(i5,' if,il=',i4,i5,' iclst:',50i4)
c      end do
      i=1
      do while (ilast(i) .gt. ifirst(i) .and. i .lt. nofcls)
        i=i+1
      end do
      nsing=0
      if (ilast(i) .eq. ifirst(i)) nsing=nofcls-i+1
      write (iout,*) 'Number of single nodes=',nsing
      nofcls_t=nofcls-nsing
      return
2000  format(' Number of nodes with fewer than',i4,' neighbors=',i5)
2001  format(' Range of the number of neighbors: [',i4,',',i4,']')
      end
