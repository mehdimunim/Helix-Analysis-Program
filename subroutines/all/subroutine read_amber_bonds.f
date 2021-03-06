      subroutine read_amber_bonds(ih,inp_top,nneig,nhneig,ineig,iatnum,
     -  read_format,lread_format,ndig,n,nr,liner,ierr,maxng,maxrec)
      dimension nneig(n),nhneig(n),ineig(maxng,n),iatnum(n)
      character*25 read_format
      character*80 liner
      dimension i20(20)
c     print *,'READ_AMBER_BONDS ih=',ih,' nr,ndig,maxng=',nr,ndig,maxng
      ierr=0
      i12=1
      i1=0
      call blankout(liner,1,80)
      do while (liner(1:5) .ne. '%FLAG')
        call blankout(liner,1,80)
        read (inp_top,1000,end=613) liner
        if (liner(1:5) .ne. '%FLAG') then
          call lastchar(liner,lc,80)
          nc=lc/ndig
          read (liner,read_format(1:lread_format),end=613)
     -      (i20(i),i=1,nc)
          do i=1,nc
            if (i12 .eq. 1) then
              i1=i20(i)/3+1
              i12=2
            else if (i12 .eq. 2) then
              i2=i20(i)/3+1
              nneig(i1)=nneig(i1)+1
              nneig(i2)=nneig(i2)+1
              ineig(nneig(i1),i1)=i2
              ineig(nneig(i2),i2)=i1
              if (iatnum(i1) .eq. 1) nhneig(i2)=nhneig(i2)+1
              if (iatnum(i2) .eq. 1) nhneig(i1)=nhneig(i1)+1
              i12=3
              if (ih .eq. 1) then
                if (iatnum(i1) .ne. 1 .and. iatnum(i2) .ne. 1) then
                  write (6,2000) i1,i2,' not'
                  ierr=1
                end if
              else 
                if (iatnum(i1) .eq. 1 .or. iatnum(i2) .eq. 1) then
                  write (6,2000) i1,i2,' '
                  ierr=1
                end if
              end if
            else if (i12 .eq. 3) then
              nr=nr+1
              call checkdim(nr,maxrec,'MAXREC',6,'number of bonds',
     -          15,0)
              i12=1
            end if
          end do
c         write (6,8767) (i,itemp1(i),itemp2(i),i=nr0+1,nr)
c8767     format(i5,' IT1,2=',2i6)
        end if
      end do
      if (i12 .ne. 1) then
        print *,'PROGRAM ERROR: bond list length is not a multiple ',
     -    'of three'
        stop
      end if
      return
613   ierr=1
      return
1000  format(a)
2000  format(' (PROGRAM?) ERROR: bond ',i7,' - ',i7,' does',a,
     -  ' involve a hydrogen')
      end
