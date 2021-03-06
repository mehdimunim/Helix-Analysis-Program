      subroutine nextnames(trajnam,ltrajnam,trajnam_n,ltrajnam_n,ntraj)
      character*200 trajnam,trajnam_n(2)
      character*80 ext
      dimension ltrajnam_n(2)
c     Find out run number/version number category and generate new name options
      do i=1,2
        trajnam_n(i)=trajnam
      end do
      lc=ltrajnam
      do while (trajnam(lc:lc) .ne. '.')
        lc=lc-1
      end do
      lext=ltrajnam-lc+1
      ext(1:lext)=trajnam(lc:ltrajnam)
      lc=lc-1
      lv2=lc
      do while (trajnam(lc:lc) .ne. '.' .and.
     -          trajnam(lc:lc) .ne. '_' .and. lc .gt. 1)
        lc=lc-1
      end do
      if (lc .ge. lv2) go to 100
      read (trajnam(lc+1:lv2),*,err=100) nr
      len_nr=lv2-lc
c     Number read
      if (trajnam(lc:lc) .eq. '.') then
c       x.nr.ext
        lr=lc+1
        ntraj=2
        lc=lv2+2
        trajnam_n(1)(lv2+1:lc)='_2'
        trajnam_n(1)(lc+1:lc+lext)=ext(1:lext)
        ltrajnam_n(1)=lc+lext
        nr=nr+1
        call writeint(trajnam_n(2),lr,nr,len)
        lc=lr-1
        trajnam_n(2)(lc+1:lc+lext)=ext(1:lext)
        ltrajnam_n(2)=lc+lext
      else
        nv=nr
c       x.$_nv.ext
        lr2=lc-1
        do while (trajnam(lc:lc) .ne. '.' .and. lc .gt. 1)
          lc=lc-1
        end do
        lr=lc+1
        if (lc .ge. lr2) go to 200
        read (trajnam(lc+1:lr2),*,err=200) nr
c       x.nr_nv.ext
        nv0=nv
        ntraj=2
        nv=nv+1
        lc=lv2-len_nr+1
        print *,'lv2,lc=',lv2,lc
        call writeint(trajnam_n(1),lc,nv,len)
        lc=lc-1
        trajnam_n(1)(lc+1:lc+lext)=ext(1:lext)
        ltrajnam_n(1)=lc+lext
        nr=nr+1
        lr0=lr
        call writeint(trajnam_n(2),lr,nr,len)
        trajnam_n(2)(lr:lr)='_'
        lr=lr+1
        call writeint(trajnam_n(2),lr,nv0,len)
        lc=lr-1
        trajnam_n(2)(lc+1:lc+lext)=ext(1:lext)
        ltrajnam_n(2)=lc+lext
      end if
      return
200   continue
c     x_nv.ext
      ntraj=2
      lc=lr2+3
      trajnam_n(1)(lr2+1:lc)='.2_'
      lc=lc+1
      call writeint(trajnam_n(1),lc,nv,len)
      lc=lc-1
      trajnam_n(1)(lc+1:lc+lext)=ext(1:lext)
      ltrajnam_n(1)=lc+lext
      lc=lv2
      nv=nv+1
      call writeint(trajnam_n(2),lc,nv,len)
      lc=lc-1
      trajnam_n(2)(lc+1:lc+lext)=ext(1:lext)
      ltrajnam_n(2)=lc+lext
      return
c     x.ext
100   ntraj=2
      lc=lv2+2
      trajnam_n(1)(lv2+1:lc)='.2'
      trajnam_n(1)(lc+1:lc+lext)=ext(1:lext)
      ltrajnam_n(1)=lc+lext
      trajnam_n(2)(lv2+1:lc)='_2'
      trajnam_n(2)(lc+1:lc+lext)=ext(1:lext)
      ltrajnam_n(2)=lc+lext
      return
      end
