      subroutine normalmodes(ncorr,inpt,nframe,inptyp,iout,iannout,
     -  ierr,index,value,ifa,ila,itemp,temp,maxt)
c     Calculate the normal modes from the residue correlation matrix
      dimension index(maxt),value(maxt),ifa(maxt),
     -  ila(maxt),itemp(maxt),temp(maxt)
      parameter (MAXPHI=400,MAX2D=5000)
      parameter (IFILL4=MAXPHI*MAXPHI*MAXPHI-
     -  (2*MAX2D*MAX2D+17*MAX2D))
      real*8 trajcorr,diag,offdiag,diagfill,drow,drowfill,cavs1,cavs2
      common /nnwork/ trajcorr(MAX2D,MAX2D),drow(MAX2D),
     -  drowfill(2,MAX2D),diag(MAX2D),offdiag(MAX2D),
     -  diagfill(MAX2D),cavs1(MAX2D),cavs2(MAX2D),
     -  row(MAX2D),fill(IFILL4)
      character*80 lineinp
c     print *,'NORMALMODES inpt=',inpt,' inptyp=',inptyp
c     Read in the covariance matrix
      ierr=1
      rewind inpt
      if (inptyp .eq. 0) then
c       Read matrix from the temp file written by Simulaid
        do jr=1,ncorr
          read (inpt,end=999) (drow(ir),ir=1,ncorr)
          do ir=1,ncorr
            trajcorr(ir,jr)=drow(ir)
          end do
        end do
      else if (inptyp .eq. 1) then
c       Binary
        read (inpt,end=100,err=100) ncorr
        do i=1,ncorr
          read (inpt,end=101,err=101) (row(j),j=1,ncorr)
          do j=1,ncorr
            trajcorr(i,j)=row(j)
          end do
        end do
      else
c       Ascii
        read (inpt,1005,end=100,err=100) lineinp
        call lastchar(lineinp,lc,80)
        ic=1
        if (lineinp(ic:ic) .eq. '#') ic=2
        read (lineinp(ic:lc),*,end=100,err=100) ncorr
        if (inptyp .eq. 2 .and. ncorr .gt. 3000) then
          print *,'Change format 2006 in subroutine normalmodes'
          stop
        end if
        do i=1,ncorr
          if (inptyp .eq. 3) then
            read (inpt,1004,end=101,err=101) (trajcorr(i,j),j=1,ncorr)
          else
            read (inpt,1006,end=101,err=101) (trajcorr(i,j),j=1,ncorr)
          end if
        end do
      end if
      ierr=0
c     write (78,*) 'ncorr=',ncorr
c     do i=1,ncorr
c       write (78,1004) (trajcorr(i,j),j=1,ncorr)
c     end do
      call dtred2(trajcorr,ncorr,MAX2D,diag,offdiag)
      call dtqli(diag,offdiag,ncorr,MAX2D,trajcorr,ierr)
      if (ierr .gt. 0) then
        write (6,1007)
        write (iout,1004)
        return
      end if
c     Print normal modes
c     Columns of trajcorr are the eigenvectors
      call askyn('Do you want to sort by eigenvalues',34,1,-1,isortev,0,
     -  0)
      call indexit(index,1,ncorr,0)
      if (isortev .eq. 1) then
        do i=1,ncorr
          value(i)=-diag(i)
        end do
        call mrgsrt(6,index,value,ncorr,ifa,ila,itemp,temp,maxt)
      end if
      if (iannout .eq. 1) then
        if (inptyp .eq. 0) then
          call write_traj_lim(iout,
     -      'Eigenvalues and eigenvectors of the covariance matrix',53,
     -      1,incr_tr,0)
          write (iout,1001) nframe
        end if
        do jjr=1,ncorr
          jr=index(jjr)
          write (iout,1002) jjr,diag(jr),(trajcorr(ir,jr),ir=1,ncorr)
        end do
      else
        do jjr=1,ncorr
          jr=index(jjr)
          write (iout,1003) diag(jr),(trajcorr(ir,jr),ir=1,ncorr)
        end do
      end if
      return
999   write (6,1000) jr
      return
100   print *,'Input covariance matrix is empty'
      return
101   print *,'Input covariance matrix read aborted at row ',i
      return
1000  format(' ERROR: residue correlation matrix of file trajcorr.mat',
     -  ' is incomplete',/,' Reading column ',i5,/,
     -  ' Normal mode calculation is skipped')
1001  format(' Number of frames used=',i8)
1002  format(i6,' Eigenvalue=',e12.5,' eigenvector:',/,(5e13.5))
1003  format(3000e13.6)
1004  format(5e13.6)
1005  format(a)
1006  format(3000e13.6)
1007  format(' Calculation aborted due to diagonalization failure')
      end
