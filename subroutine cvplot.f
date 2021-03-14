      subroutine cvplot(c,n,nslt,imaptyp,line,index,indexr,inamcol1,
     -  inamcol2,iresncol1,iresncol2,repnam,title,ltitle,ncolcode,
     -  maxcolcode,iwr,ips,iclose,maxrec,ipspage)
      dimension c(3,n),index(n),indexr(n)
      character*8 repnam
      character*(*) title
      parameter (MAXPHI=400,MAXCV=2000)
      parameter (IFILL3=MAXPHI*MAXPHI*MAXPHI-(8*MAXCV+5*MAXCV*MAXCV))
      common /nnwork/ mx(MAXCV,MAXCV),rij(3,MAXCV,MAXCV),
     -  dij(MAXCV,MAXCV),cvavfor(MAXCV),cvavback(MAXCV),cvav(MAXCV),
     -  cvrow(MAXCV),cvcol(MAXCV),ca(3,MAXCV),fill(IFILL3)
      dimension rmx(1,1),ixshuffle(MAXCV)
      real*8 rijsum(3),dijsum,forsum(MAXCV),backsum(MAXCV),dmx(1,1)
      character*4 yclab(1)
      character*47 maptyp(2)
      character*8 atnam
      character* 132 line(maxrec)
      common /graphics/ npixx,npixy,maxpixx,maxpixy,idwmain,idwplot,
     -  wx,wy,wz,wxdr
      data mxr /MAXCV/,iok /0/
      data maptyp /
     -  'Domain map based on standard circular variances',
     -  'Domain map based on weighted circular variances'/
      data nyclab /1/,lyclab /1/
      call indexit(ixshuffle,1,MAXCV,0)
      lnam=inamcol2-inamcol1+1
c     print *,'inamcol1,inamcol2,lnam=',inamcol1,inamcol2,lnam
c     print *,'renam=',repnam(1:lnam),'*'
      nres=0
      ireso=0
      icafound=0
      do ia=1,nslt
        atnam(1:lnam)=line(index(ia))(inamcol1:inamcol2)
        call leftadjustn(atnam,atnam,8)
        read (line(index(ia))(iresncol1:iresncol2),*,ERR=999) iresn
c       write (77,*) ia,iresn,' atnam=',atnam(1:lnam),'*'
        if (iresn .ne. ireso) then
c         New residue
          if (ireso .ne. 0) then
            iok=1
            if (icafound .lt. 1) then
              print *,'Residue ',ireso,' had no ',repnam(1:lnam)
              iok=0
            end if
            if (iok .eq. 1) then
              nres=nres+1
            end if
          end if
          ireso=iresn
          icafound=0
        end if
        if (atnam(1:lnam) .eq. repnam(1:lnam)) then
          if (nres .ge. mxr) then
            print *,'ERROR: program is prepared only for ',mxr,
     -        ' residues'
            stop
          end if
          icafound=1
          call trnsfr(ca(1,nres+1),c(1,ia),3)
        end if
      end do
      if (iok .eq. 1) nres=nres+1
      if (nres .eq. 0) then
        print *,'No residues containing ',repnam(1:lnam),' were found'
        stop
      end if
      print *,'Number of ',repnam(1:lnam),'-containig residues=',nres
      do ix=1,nres
        do iy=ix,nres
          d=0.0
          do k=1,3
            rij(k,iy,ix)=ca(k,iy)-ca(k,ix)
            rij(k,ix,iy)=-rij(k,iy,ix)
            d=d+rij(k,iy,ix)**2
          end do
          d=sqrt(d)
          if (imaptyp .eq. 1 .and. iy .ne. ix) then
c           Normalize the rij vectors
            do k=1,3
              rij(k,iy,ix)=rij(k,iy,ix)/d
              rij(k,ix,iy)=-rij(k,iy,ix)
            end do
            d=1.0
          end if
          dij(iy,ix)=d
          dij(ix,iy)=d
        end do
      end do
c     Internally 1-CV is stored/calculated
      cvavfor(1)=0.0
      cvavfor(nres)=0.0
      cvavback(1)=0.0
      cvavback(nres)=0.0
      do ix=1,nres
        mx(ix,ix)=9
        call zeroitd(rijsum,3)
        dijsum=0.d0
        forsum(ix)=0.d0
        do iy=ix+1,nres
          do k=1,3
            rijsum(k)=rijsum(k)+rij(k,ix,iy)
          end do
          dijsum=dijsum+dij(ix,iy)
          cv=dsqrt(rijsum(1)**2+rijsum(2)**2+rijsum(3)**2)/dijsum
          cvrow(iy)=cv
          forsum(ix)=forsum(ix)+cv
          mx(ix,iy)=float(ncolcode)*cv+1
          if (mx(ix,iy) .lt. 1) mx(ix,iy)=1
        end do
        if (nrep .le. 1 .and. iwr .gt. 0) write (iwr,2000)
     -     'For',ix,forsum(ix),(cvrow(iy),iy=ix+1,nres)
      end do
c     Now scan backwards
      do ix=nres,1,-1
        call zeroitd(rijsum,3)
        dijsum=0.d0
        backsum(ix)=0.d0
        do iy=ix-1,1,-1
          do k=1,3
            rijsum(k)=rijsum(k)+rij(k,ix,iy)
          end do
          dijsum=dijsum+dij(ix,iy)
          cv=dsqrt(rijsum(1)**2+rijsum(2)**2+rijsum(3)**2)/dijsum
          cvcol(iy)=cv
          backsum(ix)=backsum(ix)+cv
          mx(ix,iy)=float(ncolcode)*cv+1
          if (mx(ix,iy) .lt. 1) mx(ix,iy)=1
        end do
        if (nrep .le. 1 .and. iwr .gt. 0)
     -    write (iwr,2000) 'Bac',ix,backsum(ix),(cvcol(iy),iy=1,ix-1)
      end do
      do ix=2,nres-1
        cvavfor(ix)=forsum(ix)/(nres-ix)
        cvavback(ix)=backsum(ix)/(ix-1)
        cvav(ix)=(forsum(ix)+backsum(ix))/nres
      end do
      if (nrep .le. 1 .and. iwr .gt. 0)
     -  write (iwr,2001) (cvav(ix),ix=1,nres)
c      do i=1,10
c        write (6,7733) i,(mx(i,j),j=1,10)
c      end do
c7733  format(' MX i=',i3,' mx(i)=',10i3)
c     Plot map
      nrep=0
      nrep=nrep+1
      inc=max0(1,500/nres)
      ixdel=25
      iydel=115
      iytop=0
      call plotmat(ips,mx,rmx,dmx,nres,nres,0,0,0,0,1,nrep,ixdel,iydel,
     -  0,iytop,0.0,1.0,ncolcode,maxcolcode,ixdelsh,iydelsh,inc,1.0,
     -  indexr,ixshuffle,ixshuffle,title,ltitle,' ',0,0,' ',1,ca,yclab,
     -  nyclab,lyclab,mxr,1,1,MAXCV,mxr,ipspage,0)
      iydown=-80
c     Plot row averages for domain map
      if (nres .lt. 750) ixdel=ixdel+45
      scalefac=inc
      call colstrip(cvavfor,nres,ixdel,scalefac,iydown,nrep,' F',
     -  ncolcode,ips)
      call colstrip(cvavback,nres,ixdel,scalefac,iydown+18,nrep,' B',
     -  ncolcode,ips)
      call colstrip(cvav,nres,ixdel,scalefac,iydown+36,nrep,' T',
     -  ncolcode,ips)
      iydown=iydown+54
      call colcode01(ips,ixdel,iydown+5,ncolcode,nrep)
      iydown=iydown+18
c     Plot the type of the plot and close the plot/page
      if (nrep .le. 1) then
        call rgbcolor(ips,9)
        write (ips,3001) ixdel,iydown
        call psshow(ips,maptyp(imaptyp),47)
      end if
      write (ips,*) 'showpage'
      if (iclose .eq. 1) then
        close (iwr)
        close (ips)
      end if
      return
999   print *,'ERROR: illegal residue number in line',index(ia),':'
      print *,line(index(ia))(1:80)
      stop
2000  format(1x,a,i4,f10.4,(20f5.2))
2001  format(' CV column averages:',/,(8x,20f5.2))
3001  format(i4,i5,' m')
      end
