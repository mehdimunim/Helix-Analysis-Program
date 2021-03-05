      subroutine read_bondmat(rnb,inpt,lab,llab,nres,inpmat,linpmat,
     -  itemp,mxdim)
      dimension rnb(mxdim,mxdim),itemp(nres)
      character*(*) lab
      character*80 inpmat,prompt,line
      do ir=1,nres
        do jr=1,nres
          rnb(ir,jr)=0.0
        end do
      end do
      prompt(1:llab)=lab(1:llab)
      prompt(llab+1:llab+28)=' residue-residue bond matrix'
      lprompt=llab+28
      linpmat=0
      call openfile(inpt,0,prompt,lprompt,'old',inpmat,linpmat,notfnd,
     -  0,1,1,0,0)
      call indexit(itemp,1,nres,0)
c     Read pointer if necessary
      do while (line(1:19) .ne. ' Average number of ' .and.
     -  line(1:12) .ne. ' iy(orig):')
        call blankout(line,1,80)
        read (inpt,1000,end=999) line
      end do
      ixf=0
      ixl=1
      if (line(1:12) .eq. ' iy(orig):') then
c       Read indices
        ifr=1
        ilr=10
        call lastchar(line,lc,80)
        if (lc .lt. 73) ilr=ilr-(73-lc)/5
        read (line(24:73),1002,end=888,err=888) (itemp(i),i=ifr,ilr)
        do while (lc .eq. 73)
          ifr=ilr+1
          ilr=ilr+10
          call blankout(line,1,80)
          read (inpt,1000,end=999) line
          call lastchar(line,lc,80)
          if (lc .lt. 73) ilr=ilr-(73-lc)/5
          read (line(24:73),1002,end=888,err=888) (itemp(i),i=ifr,ilr)
        end do
        do while (line(1:19) .ne. ' Average number of ')
          read (inpt,1000,end=999) line
        end do
      end if
      iymin=99999
      iymax=0
      read (inpt,1000,end=999) line
      do while (.true.)
        ixf=1
        ixl=0
        if (line(10:18) .eq. 'iy(orig)=') then
c         Start of one column
          read (line(19:23),*,end=888,err=888) iyorig
          read (line(5:8),*,end=888,err=888) iy
          if (itemp(iy) .ne. iyorig) then
            write (6,1003) iy,iyorig,itemp(iy)
            stop
          end if
          if (iy .lt. iymin) iymin=iy
          if (iy .gt. iymax) iymax=iy
          ixf=1
          ixl=ixf+9
          read (line(25:74),1001,end=888,err=888)
     -      (rnb(itemp(irx),iyorig),irx=ixf,ixl)
c           write (88,7211) iy,iyorig,ixf,ixl,
c    -        (rnb(itemp(irx),iyorig),irx=ixf,ixl)
c7211       format(' iy,iyorig=',2i5,' ixf,ixl=',2i5,' r=',10f7.2)
          read (inpt,1000,end=999) line
          do while(line(10:18) .ne. 'iy(orig)=' .and.
     -             line(1:17) .ne. ' Fraction of time')
            ixf=ixl+1
            ixl=ixf+9
            read (line(25:74),1001,end=888,err=888)
     -        (rnb(itemp(irx),iyorig),irx=ixf,ixl)
c           write (88,7212) ixf,ixl,
c    -        (rnb(itemp(irx),iyorig),irx=ixf,ixl)
c7212       format(' ixf,ixl=',2i5,' r=',10f7.2)
            read (inpt,1000,end=999) line
          end do
          if (line(1:17) .eq. ' Fraction of time') go to 999
        end if
      end do
999   close (inpt)
      write (6,1004) iymin,iymax
      return
888   print *,'ERROR reading record ixf,ixl=',ixf,ixl
      print *,line(1:66)
      close (inpt)
      return
1000  format(a)
1001  format(10f5.2)
1002  format(10i5)
1003  format(' PROGRAM ERROR: iy=',i5,' iyorig=',i5,' iyorig(iy)=',i5)
1004  format(' Bond matrix info read for indices [',i5,',',i5,']')
      end
