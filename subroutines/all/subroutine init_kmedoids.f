      subroutine init_kmedoids(nofcls,index2d,ifirst,ilast,icent,rmsd2d,
     -  n0,n,nnode,ierr,iverb,mx2d)
      dimension index2d(mx2d),icent(mx2d),ifirst(mx2d),ilast(mx2d),
     -  rmsd2d(mx2d,mx2d)
      integer*4 ixo
      common /rangen/ ixo
      dimension ran(1)
      character*31 qcent
c     print *,'INIT_KMEDOIDS n0,n,mx2d=',n0,n,mx2d
      if (nofcls .le. 0 .or. nofcls .gt. n-n0+1) then
        print *,'PROGRAM ERROR: invalid number of clusters requested:',
     -    nofcls
        nofcls=1
        ifirst(1)=n0
        ilast(1)=n
        ierr=1
        return
      end if
      nnode=n-n0+1
      if (nnode .lt. nofcls) then
        print *,'ERROR: can not generate ',nofcls,' clusters from',nnode
        ierr=1
        return
      end if
      ierr=0
      call askyn('Do you want to specify initial cluster centers',46,
     -  1,-1,ireadct,92,0)
      if (ireadct .eq. 1) then
        do k=1,nofcls
          write (qcent,1001) k
          call getint(qcent,31,999999,1,n,icent(k),00)
        end do
      else
        if (ixo .ne. 1237) then
          call askyn(
     -      'Do you want to reinitialize the random-number seed',
     -      50,1,1,ireinit,92,0)
          if (ireinit .eq. 1) call randpx_init(1357)
        end if
        call randpx(1,ran)
        icent(1)=n0+ran(1)*nnode
        ix=0
        do k=2,nofcls
          rminmax=0.0
          do i=n0,n
            ii=index2d(i)
            rmin=10000.0
            do j=1,k-1
              if (i .eq. icent(j)) then
                rmin=0.0
              else
                if (rmin .gt. rmsd2d(ii,index2d(icent(j))))
     -            rmin=rmsd2d(ii,index2d(icent(j)))
              end if
            end do
c           Select cent(k) to have the largest rmin
            if (rmin .gt. rminmax) then
              icent(k)=i
              rminmax=rmin
            end if
          end do
        end do
      end if
      if (iverb .gt. 0) write (6,1000) (icent(k),k=1,nofcls)
      return
1000  format(' Initial centers chosen=',(10i5))
1001  format(' Initial center for cluster #',i3)
      end
