      subroutine findtcent(ineig,nneig,iparent,list,n0,n,icent,
     -  iout,maxneig,maxat)
c#    MMC routine 463 lstmod: 01/20/05
c*****Find the topology center of a molecule
      dimension ineig(maxneig,maxat),nneig(maxat),list(maxat),
     -  iparent(maxat)
      data ilist /0/
c     print *,'FINDTCENT n0,n,maxneig,maxat=',n0,n,maxneig,maxat
      if (n .eq. n0) then
        icent=n0
        return
      end if
      do ir=1,2
        nsteps=0
        do i=n0,n
          iparent(i)=0
        end do
        if (ir .eq. 1) then
c         First grow neighbour list from n0
          ic=n0
        else
c         Now grow from list(ilist)
          ic=list(ilist)
        end if
        iparent(ic)=ic
        ilistf=1
        ilistl=1
        list(ilistf)=ic
        do while (ilistl .ge. ilistf)
          nsteps=nsteps+1
          ilist=ilistl
          do il=ilistf,ilistl
            ic=list(il)
            do in=1,nneig(ic)
              ia=ineig(in,ic)
              if (ia .ge. n0 .and. ia .le. n) then
c               Consider only atoms within the range
                if (iparent(ia) .eq. 0) then
                  iparent(ia)=ic
                  ilist =ilist+1
                  list(ilist)=ia
                end if
              end if
            end do
          end do
          ilistf=ilistl+1
          ilistl=ilist
        end do
c       Now list(ilist) is one of the farthest from ic
c       write (77,1000) n0,n,ir,list(ilist),nsteps
c1000  format(' FINDCENT: Center search for atoms ',i5,' - ',i5,
c     -  ': run',i2,' End point: ',i5,' Steps: ',i4)
      end do
c     nsteps is the number of vertices on the longest path, backtrack
c     nsteps/2 to get the center
      ic=list(ilist)
      do ia=1,nsteps/2
        ic=iparent(ic)
      end do
      icent=ic
      if (iout .gt. 0) write (iout,1001) icent,n0,n
      return
1001  format(' Center found:',i5,' in atom range [',i5,',',i6,']')
      end
