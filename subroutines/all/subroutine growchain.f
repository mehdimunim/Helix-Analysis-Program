      subroutine growchain(n0,n,n00,iaf,ial,nsteps,nneig,ineig,
     -  list,iparent,ichain,imax,maxng)
      dimension nneig(n),ineig(maxng,n),list(n),iparent(n),ichain(n)
c     write (77, *) 'GROWCHAIN n0,n,n00,iaf=',n0,n,n00,iaf
c     write (77, *) 'GROWCHAIN nforbid,iaforbid=',nforbid,iaforbid
      nsteps=0
      imax=0
      call zeroiti(iparent,n0-1,n)
      if (n00 .ne. 0) iparent(n00)=n00
      iparent(iaf)=iaf
      ic=iaf
      ilistf=1
      ilistl=1
      list(ilistf)=ic
      imax=ic
      do while (ilistl .ge. ilistf)
        nsteps=nsteps+1
        ilist=ilistl
        do il=ilistf,ilistl
          ic=list(il)
          do in=1,nneig(ic)
            ia=ineig(in,ic)
            if (ia .ge. n0 .and. ia .le. n) then
c             Consider only atoms within the range
              if (iparent(ia) .eq. 0) then
                iparent(ia)=ic
                ilist =ilist+1
                list(ilist)=ia
                if (ia .gt. imax) imax=ia
              end if
            end if
          end do
        end do
        ilistf=ilistl+1
        ilistl=ilist
      end do
c     nsteps is the number of vertices on the longest path, backtrack
      ial=list(ilist)
      ic=ial
      ichain(1)=ic
      do ia=2,nsteps+1
        ic=iparent(ic)
        ichain(ia)=ic
      end do
c      write (77,1000) n00,iaf,(ichain(ia),ia=1,nsteps+1)
c1000  format(' n00,iaf=',2i4,(' chain:',25i4))
      return
      end
