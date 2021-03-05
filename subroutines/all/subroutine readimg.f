      subroutine readimg(cell,ncell)
      dimension cell(3,27)
      character*1 xyz
      common /axislab/ xyz(3)
      character*1 abc,digits,hexdigits
      common /charactersets/ ihex(25),abc(62),digits(14),hexdigits(25)
      character*1 ans
      character*80 line
      character*200 imgfile
      dimension scale(3)
c     print *,'READIMG cell=',scale
      call quiz(ans,ians,' ','image',5,'PBC cell type',13,0,5,6,0)
      lenfilename=0
      call openfile(30,0,'image',5,'old',imgfile,lenfilename,notfnd,
     -  1,1,1,0,0)
      if (ans .eq. 'c') then
c       Charmm
        do k=1,3
          call getreal(
     -    'Charmm '//xyz(k)//' axis image parameter (@'//abc(58+k)//')',
     -      34,999999.0,scale(k),1,0)
        end do
        ncell=1
        call zeroit(cell,3)
        rewind 30
        line(1:3)='   '
        do while (line(1:3) .ne. 'END' .and. line(1:3) .ne. 'end')
          read (30,888,end=999) line
c         write (6,888) line
          if (line(1:3) .eq. 'END') then
            print *,'Read ',ncell-1,' images - (0,0,0) was added'
c           write (6,1003) (ic,(cell(k,ic),k=1,3),ic=1,ncell)
          else
            if (line(1:5) .eq. 'TRANS') then
              ncell=ncell+1
              read (line(7:36),900,err=998) (cell(k,ncell),k=1,3)
              do k=1,3
                cell(k,ncell)=-scale(k)*cell(k,ncell)
              end do
            end if
          end if
        end do
      else
c       Simulaid format
        do k=1,3
          call getreal('Multiplying factor for the '//xyz(k)//' axis',
     -      33,1.0,scale(k),1,37)
        end do
        ncell=0
        do while (.true.)
          read (30,1001,err=998,end=100) (cell(k,ncell+1),k=1,3)
          ncell=ncell+1
          do k=1,3
            cell(k,ncell)=cell(k,ncell)*scale(k)
          end do
        end do
100     print *,'Read ',ncell,' cell centers - (0,0,0) is assumed ',
     -    'to be listed'
      end if
      return
998   write (6,1000) ncell,line(7:36)
      stop
999   print *,'ERROR: Image file is not terminated properly'
      stop
900   format(3f10.0)
888   format(a80)
c1003  format(' The PBC cell centers:',/,30(i5,3f10.4,/))
1000  format(' ERROR: invalid cell coordinate input for cell ',i3,
     -  ':',/,8x,a)
1001  format(3f10.5)
      end
