      subroutine read_distmat(iw0,n,rijrmin,rijrmax,ix,label,inpfile,
     -  iuout,mx2d)
      dimension ix(mx2d)
      character*500 label(mx2d)
      character*500 line
      character*(*) inpfile
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (4*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
     -  ihbpair(2,MAXBONDS),a1(MAXBONDS),a2(MAXBONDS),fill(IFILL2)
      dimension n10(10)
      llen=500
      namleni=0
      call openfile(iw0,0,'distance matrix file',20,'old',inpfile,
     -  namleni,notfnd,0,1,1,0,0)
      call askyn('Do you want to limit the label to its last word',47,
     -  1,-1,lastword,0,0)
      read (iw0,*,err=666) n
      call checkdim(n,MAX2D,'MAX2D',5,'matrix dimension',16,0)
      if (n .lt. 1) then
        write (6,1005) n
        stop
      else
        write (6,1003) n
        write (iuout,1003) n,' ',inpfile(1:namleni)
      end if
      do i=1,n
        call blankout(line,1,llen)
        read (iw0,2000,end=999) line
        call lastchar(line,lline,llen)
        if (lline .eq. llen) write (6,1004) llen,i,line
        ic=1
        if (line(1:1) .eq. ' ') call nextchar(line,ic,llen)
        ic1=ic
        call nextblank(line,ic,llen)
c       print *,'ic1,ic=',ic1,ic
c       print *,'line=',line(1:30)
        read (line(ic1:ic-1),*,err=777) ix(i)
        call blankout(label(i),1,llen)
        if (lastword .eq. 0) then
          call nextchar(line,ic,llen)
          label(i)(1:lline-ic+1)=line(ic:lline)
        else
          call laststring(line,icf,icl,lline,llen)
          label(i)(1:icl-icf+1)=line(icf:icl)
        end if
c       print *,'i=',i,' ix=',ix(i),' label=',label(i)(1:lline-ic+1)
        read (iw0,*,end=999,err=888) (rmsd2d(i,j),j=1,n)
      end do
      close (iw0)
c     do i=1,n
c       print *, ix(i),label(i)
c       print *, (rmsd2d(i,j),j=1,n)
c     end do
      rijrmin=1000000.0
      rijrmax=-rijrmin
      nerr=0
      do i=1,n
        do j=1,i-1
          if (rmsd2d(i,j) .lt. rijrmin) rijrmin=rmsd2d(i,j)
          if (rmsd2d(i,j) .gt. rijrmax) rijrmax=rmsd2d(i,j)
          if (rmsd2d(i,j) .ne. rmsd2d(j,i)) then
            nerr=nerr+1
            if (nerr .le. 10)
     -         write (6,1001) i,j,rmsd2d(i,j),j,i,rmsd2d(j,i)
          end if
        end do
      end do
      write (6,1002) rijrmin,rijrmax
      call transform_mat(rmsd2d,n,mx2d,rijrmin,iuout)
      call transform_dist(rijrmin,rijmin)
      call transform_dist(rijrmax,rijmax)
      if (rijmin .gt. rijmax) then
        x=rijmax
        rijmax=rijmin
        rijmin=x
      end if
      write (6,1006) rijmin,rijmax
      call zeroiti(n10,0,10)
      rijrange=rijmax-rijmin
      nz=0
      do i=1,n
        do j=1,i-1
          ixd=10.0*(rmsd2d(i,j)-rijmin)/rijrange
          if (ixd .lt. 10) ixd=ixd+1
          n10(ixd)=n10(ixd)+1
          if (rmsd2d(i,j) .eq. 0.0) nz=nz+1
        end do
      end do
      write (6,1007) n10
      if (nz .gt. 0) then
        write (6,1008) nz
        call askyn('Do you want a list of identical pairs',37,
     -    1,-1,ilistdup,0,0)
        if (ilistdup .gt. 0) then
          do i=1,n
            do j=1,i-1
              if (rmsd2d(i,j) .eq. 0.0) then
                call lastchar(label(i),lci,500)
                call lastchar(label(j),lcj,500)
                if (lci .eq. lcj) then
                  if (label(i)(1:lci) .eq. label(j)(1:lcj)) then
                    write (6,1009) i,label(i)(1:lci),j,label(j)(1:lcj)
                  end if
                end if
              end if
            end do
          end do
        end if
      end if
      return
666   print *,' ERROR: invalid matrix dimension:'
      rewind iw0
      read (iw0,2000,end=999) line
      print *,line(1:20)
      stop
777   print *,' ERROR: invalid sequence number for item ',i,':',
     -  line(ic1:ic-1)
      stop
888   print *,' ERROR: invalid number read for item ',i
      stop
999   print *,' ERROR: matrix input ended unexpectedly'
      stop
1001  format(' Symmetry ERROR: R[',i5,',',i5,']=',f10.5,' and ',
     -  'R[',i5,',',i5,']=',f10.5)
1002  format(' Range of matrix values: [',f10.5,',',f10.5,']',/,
     -  ' NOTE: if the input matrix entries represent similarities',/,
     -  '       then you have to use one of the (1-r(i,j)) transform',
     -  ' options')
1003  format(' Distance matrix dimension=',i5,a,/,' Read from file ',a)
1004  format(' WARNING: label is longer than ',i4,' characters in row',
     -  i6,':',/,1x,a)
1005  format(' ERROR: Input matrix dimension is non=positive:',i9)
1006  format(' Transformed distance range: [',f10.5,',',f10.5,']')
1007  format(' Distribution of the distances:',/,10i7)
1008  format(' NOTE: there are ',i6,' identical item pairs')
1009  format(' #',i5,' (',a,') is identical with #',i5,' (',a,')')
2000  format(a)
      end
