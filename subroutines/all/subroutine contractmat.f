      subroutine contractmat(rij,nxo,nyo,nxn,nyn,navg,ixx,ixy,
     -  itemp1,itemp2,temp,isort,iouttemp,ndim)
      dimension rij(ndim,ndim),ixx(ndim),ixy(ndim),
     -  itemp1(ndim),itemp2(ndim),temp(ndim)
      character*9 filename
      navgdef=max0(1,max0(nxo,nyo)/1000)
c     print *,'CONTRACTMAT nxo,nyo,isort=',nxo,nyo,ndim,isort
      call getint(
     -  'Number of x and y entries to average in the 2D-RMSD plot',56,
     -  navgdef,1,min0(nxo,nyo),navg,122)
c     print *,'NAVG=',navg
      if (navg .eq. 1) then
        write (6,*) 'All matrix elements will be plotted'
        nxn=nxo
        nyn=nyo
        return
      else
        ntry=0
        iopen=1
        filename='rij__.tmp'
        do while (ntry .lt. 10 .and. iopen .gt. 0)
          write (filename(5:5),1000) ntry
          open(unit=iouttemp,status='new',file=filename,iostat=iopen,
     -      form='unformatted')
          ntry=ntry+1
        end do
        if (iopen .gt. 0) then
          print *,'PROGRAM ERROR: could not open rij.tmp'
          stop
        end if
        write (iouttemp) ((rij(i,j),i=1,nxo),j=1,nyo)
        if (isort .gt. 0) then
c         Sort columns first
          do iy=1,nyo
            do i=1,nxo
              temp(i)=rij(ixx(i),iy)
            end do
            call trnsfr(rij(1,iy),temp,nxo)
          end do
          call trnsfi(itemp1,ixy,nyo)
          do i=1,nyo
            itemp2(itemp1(i))=i
          end do
          do iy=1,nyo
            call trnsfr(temp,rij(1,itemp1(iy)),nxo)
            call trnsfr(rij(1,itemp1(iy)),rij(1,iy),nxo)
            call trnsfr(rij(1,iy),temp,nxo)
            itemp1(itemp2(iy))=itemp1(iy)
            itemp2(itemp1(itemp1(iy)))=itemp1(iy)
          end do
        end if
        nxn=(nxo-1)/navg+1
        nyn=(nyo-1)/navg+1
        do ix=1,nxn
          do iy=1,nyn
            sum=0.0
            ixf0=(ix-1)*navg
            iyf0=(iy-1)*navg
            ixl=min0(nxo,ixf0+navg)
            iyl=min0(nyo,iyf0+navg)
            do jx=ixf0+1,ixl
              do jy=iyf0+1,iyl
                sum=sum+rij(jx,jy)
              end do
            end do
            if (ixl .eq. ixf0 .or. iyl .eq. iyf0) then
              write (6,9781) ixl,ixf0,iyl,iyf0
9781          format(' ixl,ixf0=',2i6,' iyl,iyf0=',2i6)
            end if
            rij(ix,iy)=sum/((ixl-ixf0)*(iyl-iyf0))
          end do
        end do
      end if
      return
1000  format(i1)
      end
