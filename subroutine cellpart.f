      subroutine cellpart(c,ian,itemp1,na1,na,padding,spacing,corner,
     -  ecell,edge,xyzmin,xyzmax,vtot,nxyz,ixyz,icellno,itemp2,index,
     -  ifirst,ilast,ntotcell,levout,iout,maxrec)
      dimension c(3,na),ian(na),itemp1(maxrec),corner(3),cent(3),
     -  edge(3),ecell(3),xyzmin(3),xyzmax(3),nxyz(3),ixyz(3),
     -  icellno(maxrec),itemp2(maxrec),index(maxrec),ifirst(maxrec),
     -  ilast(maxrec)
c     Partition the solute atoms into cells for fast distance calculation
c     print *,'CELLPART NA1,NA,SPACING,PADDING=',na1,na,spacing,padding
      call extension(c,ian,0,na1,na,xyzmin,xyzmax,cent,0,0,v)
      do k=1,3
        edge(k)=xyzmax(k)-xyzmin(k)+2.0*padding
        corner(k)=xyzmin(k)-padding
        nxyz(k)=int(edge(k)/spacing)
        if (nxyz(k) .eq. 0 .or. edge(k)-nxyz(k)*ecell(k) .ge. 0.5)
     -    nxyz(k)=nxyz(k)+1
        ecell(k)=(edge(k)+0.0001)/nxyz(k)
      end do
      ntotcell=nxyz(1)*nxyz(2)*nxyz(3)
      if (ntotcell .gt. maxrec) then
        write (6,7000) ntotcell,maxrec,ntotcell
        stop
      end if
      vtot=edge(1)*edge(2)*edge(3)
c     Establish the cell of each solute atom
      nha=0
      do ia=na1,na
        if (ian(ia) .gt. 1) then
          do k=1,3
            ixyz(k)=(c(k,ia)-corner(k))/ecell(k)
          end do
          nha=nha+1
          icellno(nha)=1+ixyz(1)+nxyz(1)*ixyz(2)+nxyz(1)*nxyz(2)*ixyz(3)
          index(nha)=ia
          if (levout .gt. 3) write (iout,7001) ia,nha,icellno(nha),ixyz
        end if
      end do
      if (levout .gt. 2) write (iout,*) 'na,nha=',na,nha
      call mrgsrti(6,index,icellno,nha,ifirst,ilast,itemp1,itemp2,
     -  maxrec)
      if (levout .gt. 3)
     -  write (iout,7002) (ia,icellno(ia),index(ia),ia=1,nha)
      call zeroiti(ifirst,0,ntotcell)
c     Establish the limits of cells
      ia=1
      icell=icellno(ia)
      ifirst(icell)=ia
      do while (ia .lt. nha)
        ia=ia+1
        if (icellno(ia) .gt. icell) then
          ilast(icell)=ia-1
          icell=icellno(ia)
          ifirst(icell)=ia
        end if
      end do
      ilast(icell)=nha
c     For empty cells, ifirst(icell)=0!
c     print *,'return CELLPART NA=',na
      return
7000  format(' ERROR: the number of solute cells (',i8,') exceeds the',
     -  ' maximum number',/,' of records (',i8,')',/,' Recompile with ',
     -  'the parameter MAXREC > ',i8,' or increase the grid spacing')
7001  format(i5,' nha=',i4, ' icellno=',f8.1,' ixyz=',3i6)
7002  format(i5,' icellno=',f8.1,' index=',i6)
      end
