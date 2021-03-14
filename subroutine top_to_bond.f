      subroutine top_to_bond(nntyp,nneig,nhneig,ineig,iatnum,n,inp,
     -  itemp1,itemp2,maxng,maxrec)
      dimension nneig(n),nhneig(n),ineig(maxng,n),iatnum(n),
     -  itemp1(maxrec),itemp2(maxrec)
      dimension i8(8)
      character*1 ansrun
      character*25 read_format
      character*80 inpfiletmp,liner
c     print *,'TOP_TO_BOND n=',n,' maxng=',maxng
      if (nntyp .eq. 0) then
        call quiz(ansrun,nntyp,' ',' ',0,'Bond information source',23,
     -    0,5,6,00)
        if (nntyp .eq. 1) return 
      end if
      if (nntyp .eq. 2) then
c       Get nnlist from a Charmm PSF file
        if (inp .eq. 0) then
          inp_top=62
          nerr=1
          do while (nerr .gt. 0)
            namlenn=0
            call openfile(inp_top,0,'PSF',3,'old',inpfiletmp,namlenn,
     -        notfnd,3,1,1,0,0)
            call find_n_psf(inp_top,liner,nerr,n,natspsf,'!NATOM',6)
            if (n .ne. natspsf) then
              print *,'Number of atoms in the PSF file (',natspsf,
     -          ' differs from number of atoms (',n,')'
              if (natspsf .gt. n) then
                call askstop(1)
              else
                stop
              end if
            end if
          end do
        else
          inp_top=inp
          rewind inp_top
        end if 
        call find_n_psf(inp_top,liner,nerr,0,natbond,'!NBOND',6)
        call checkdim(natbond,maxrec,'MAXREC',6,'number of bonds',15,0)
        call zeroiti(nneig,0,natspsf)
        call zeroiti(nhneig,0,natspsf)
        nwr=natbond
        do while (nwr .gt. 0)
          nb=min0(4,nwr)
          read (inp_top,1108,err=605) (i8(i),i=1,2*nb)
605       do i=1,nb
            i1=i8(2*(i-1)+1)
            i2=i8(2*i)
c           print *,'BOND:',i1,' - ',i2
            nneig(i1)=nneig(i1)+1
            nneig(i2)=nneig(i2)+1
            ineig(nneig(i1),i1)=i2
            ineig(nneig(i2),i2)=i1
            if (iatnum(i1) .eq. 1) nhneig(i2)=nhneig(i2)+1
            if (iatnum(i2) .eq. 1) nhneig(i1)=nhneig(i1)+1
          end do
          nwr=nwr-nb
        end do
        print *,'Number of bonds read:',natbond
        do i=1,n
c         write (06,9782) i,iatnum(i),(ineig(j,i),j=1,nneig(i))
c9782     format(i4,' iatno=',i2,' in=',10i4)
        end do
      else if (nntyp .eq. 3) then
c       Get nnlist from an Amber .top file
611     if (inp .eq. 0) then
          inp_top=62
          notfnd=1
          do while (notfnd .gt. 0)
            namlenn=0
            call openfile(inp_top,0,'TOP',3,'old',inpfiletmp,namlenn,
     -        notfnd,3,1,1,0,0)
          end do
        else
          inp_top=inp
          rewind inp_top
        end if
        call find_ambertyp(inp_top,'%FLAG BONDS_INC_HYDROGEN',24,
     -    read_format,lread_format)
        call zeroiti(nneig,0,n)
        call zeroiti(nhneig,0,n)
        ierrtop=1
        ic=2
        do while (idigit(read_format(ic:ic),1) .eq. 1)
          ic=ic+1
        end do
        read (read_format(2:ic-1),*,err=613,end=613) ncol
        ic=ic+1
        icc=ic
        do while (idigit(read_format(ic:ic),1) .eq. 1)
          ic=ic+1
        end do
        read (read_format(icc:ic-1),*,err=613,end=613) ndig
        nr=0
        call read_amber_bonds(1,inp_top,nneig,nhneig,ineig,iatnum,
     -    read_format,lread_format,ndig,n,nr,liner,ierr,maxng,maxrec)
        if (liner(1:28) .eq. '%FLAG BONDS_WITHOUT_HYDROGEN') then
          read (inp_top,1000,end=613)
        else
          call find_ambertyp(inp_top,'%FLAG BONDS_WITHOUT_HYDROGEN',28,
     -      read_format,lread_format)
        end if
        call read_amber_bonds(2,inp_top,nneig,nhneig,ineig,iatnum,
     -    read_format,lread_format,ndig,n,nr,liner,ierr,maxng,maxrec)
        natbond=nr
        print *,'Number of bonds read:',natbond
        ierrtop=0
613     if (ierrtop .eq. 1) then
          print *,'TOP file had some error'
          go to 611
        end if
      else
        print *,'PROGRAM ERROR: invalid NN list code:',nntyp
        stop
      end if
      if (inp .eq. 0) close (inp_top)
      return
1000  format(a)
1108  format(8i10)
      end
