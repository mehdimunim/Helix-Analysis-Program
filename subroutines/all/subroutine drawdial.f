      subroutine drawdial(iout,edge,n,incr,idial,lab,llab,nrdiv,ndprow,
     -  ifirst,ilast,incrid,iconndial,ioutpr,mappdf,ipdfgrd,pi)
      character*(*) lab
      parameter (MAXFRAMES=50000,MAXCOPY=600)
      common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
     -  increment,increment2,res(2,MAXFRAMES,MAXCOPY),
     -  x0(MAXCOPY),y0(MAXCOPY),nxselres,ixselres(MAXCOPY)
      real*8 ddd,cosav,sinav
      dimension npdf(360),pdf(360)
      data riprev /0.0/,angprev /0.0/
c     Draw a 'dial' into the square of edge edge, with upper left corner at
c     (x0,y0), title lab(1:llab)
c     print *,'DRAWDIAL n=',n,' incrid,ll=',incrid,llab,' idial=',idial
c     print *,'DRAWDIAL ioutpr,mappdf=',ioutpr,mappdf
      pi2=2.0*pi
      rdtodg=180.0/pi
      lyshift= 95
      cx=x0(idial)+edge*0.50
c     if (ndprow .gt. 4) ldown=8
      cy=y0(idial)-edge*0.50+lyshift
      if (mappdf .eq. 0) then
        r=0.45*edge
        r0=0.04*edge
        r1=0.41*edge
      else
        r=0.40*edge
        r0=0.03*edge
        r1=0.37*edge
      end if
      write (iout,*) 'np'
c     Print 0,90,180,270 marks
      if (ndprow .lt. 3) then
        lfont=14
        lfshift=5
      else if (ndprow .lt. 5) then
        lfont=10
        lfshift=0
      else
        lfont=8
        lfshift=0
      end if
      write (iout,1005) 'Symbol',lfont
      if (ndprow .lt. 4 .or. ilast .eq. 1) then
        write (iout,1000) cx+r+2.0,cy-2.0
        call psshow(iout,'0',1)
      end if
      write (iout,1000) cx-5.0,cy+r+3.0+lfshift
      call psshow(iout,'p/2',3)
      if (ndprow .lt. 4 .or. ifirst .eq. 1) then
        write (iout,1000) cx-r-7.0-lfshift,cy-2.0
        call psshow(iout,'p',1)
      end if
      write (iout,1000) cx-8.0,cy-r-10.0-lfshift
      call psshow(iout,'3p/2',4)
      lfont=12
      if (ndprow .gt. 4) lfont=8
      write (iout,1005) 'Helvetica',lfont
      write (iout,1000) x0(idial)+5.0,y0(idial)+lyshift
      call psshow(iout,lab,llab)
      write (iout,*) 'sk'
c     Draw the circle
      write (iout,*) 'np'
      write (iout,1000) cx+r,cy
      write (iout,1003) cx,cy,r,0.0,360.0
      if (nrdiv .gt. 1) then
        write (iout,*) 'sk'
        write (iout,*) '[2] 0 setdash'
        do ir=1,nrdiv-1
          write (iout,1003) cx,cy,r0+r1*float(ir)/float(nrdiv),
     -      0.0,360.0
        end do
      end if
      write (iout,*) 'sk'
      write (iout,*) '[] 0 setdash'
      write (iout,*) 'np'
c     Draw ticks
      do ia=1,36
        ang=float(ia*10)/rdtodg
        rx=sin(ang)*r
        ry=cos(ang)*r
        write (iout,1000) cx+rx,cy+ry
        if (mod(ia,3) .eq. 0) then
          write (iout,1001) cx+0.94*rx,cy+0.94*ry
        else
          write (iout,1001) cx+0.97*rx,cy+0.97*ry
        end if
      end do
c     Draw axes
      write (iout,1000) cx-r,cy
      write (iout,1001) cx+r,cy
      write (iout,*) 'sk'
      write (iout,*) 'np'
      write (iout,1000) cx,cy+r
      write (iout,1001) cx,cy-r
      write (iout,*) 'sk'
c     Draw inner disk and initial postion line
c     write (iout,*) '0.7 0.7 0.7 setrgbcolor'
      write (iout,*) '0.6 0.6 0.6 setrgbcolor'
      write (iout,*) 'np'
      write (iout,1003) cx,cy,r0,0.0,360.0
      write (iout,*) 'fill'
      write (iout,*) 'sk'
      call rgbcolor(iout,-4)
      write (iout,1004) 2
      write (iout,*) 'np'
      write (iout,1000) cx,cy
      write (iout,1001) cx+r0*res(1,1,incrid+idial),
     -  cy+r0*res(2,1,incrid+idial)
      write (iout,*) 'sk'
      anggrid=5.0
      cosav=0.d0
      sinav=0.d0
      angmin=pi2
      angmax=0.0
      incrloop=max0(1,incr)
c     write (77,*) 'DRAWDIAL: ',lab(1:llab)
      rdot=0.5
      if (n .lt. 100) rdot=1.0
      if (n .lt. 10) rdot=1.5
      call rgbcolor(iout,-6)
      do i=1,n,incrloop
        if (incrloop .eq. 1) then
          ang=dacoscheck(ddd,res(1,i,incrid+idial),0,ioutpr,'DRAWDIAL')
          if (res(2,i,incrid+idial) .lt. 0.0) ang=pi2-ang
          if (ang. gt. angmax) angmax=ang
          if (ang. lt. angmin) angmin=ang
        else
c         First average over incrloop steps
          cosavinc=0.0
          sinavinc=0.0
          do j=i,i+min0(n,incrloop-1)
            cosavinc=cosavinc+res(1,j,incrid+idial)
            sinavinc=sinavinc+res(2,j,incrid+idial)
            ang=dacoscheck(ddd,res(1,i,incrid+idial),0,ioutpr,
     -        'DRAWDIAL')
            if (res(2,i,incrid+idial) .lt. 0.0) ang=pi2-ang
            if (ang. gt. angmax) angmax=ang
            if (ang. lt. angmin) angmin=ang
          end do
          cosavinc=cosavinc/incrloop
          sinavinc=sinavinc/incrloop
          sqsum=sqrt(cosavinc**2+sinavinc**2)
          cosavinc=cosavinc/sqsum
          sinavinc=sinavinc/sqsum
          ang=dacoscheck(ddd,cosavinc,0,iout,'DRAWDIAL')
          if (sinavinc .lt. 0.0) ang=-ang
        end if
        ri=r0+r1*float(i)/float(n)
        if (iconndial .eq. 1) then
          if (i .eq. 1) then
            write (iout,1004) 1
            write (iout,*) 'np'
            write (iout,1000) cx+cos(ang)*r0,cy+sin(ang)*r0
          else
            angdev=(ang-angprev)
            if (angdev .gt. 0.0) then
              angdevrev=angdev-2.0*pi
            else
              angdevrev=angdev+2.0*pi
            end if
            if (abs(angdevrev) .lt. abs(angdev)) angdev=angdevrev
c           if (angdev .gt. pi) angdev=angdev-2.0*pi
c           if (angdev .lt. -pi) angdev=angdev+2.0*pi
            if (abs(angdev)*ri .gt. anggrid) then
              nd=abs(angdev)*ri/anggrid
c           write (77,8788) ang*rdtodg,angprev*rdtodg,angdev*rdtodg,
c    -        angdevrev*rdtodg,nd
c8788       format(' ang,prev=',2f7.1,' angdev,rev=',2f7.1,' nd=',i4)
              if (i .gt. n/2) nd=2*nd
c             if (lab(1:1) .eq. 'P') print *,'i,ri,angdev,nd=',
c    -            i,ri,angdev,nd
              do id=1,nd
                angi=angprev+angdev*float(id)/float(nd)
                rii=riprev+(ri-riprev)*float(id)/float(nd)
                write (iout,1001) cx+cos(angi)*rii,cy+sin(angi)*rii
              end do
            else
              write (iout,1001) cx+cos(ang)*ri,cy+sin(ang)*ri
            end if
          end if
        else
c         Just draw dots
          ri=r0+r1*float(i)/float(n)
          write (iout,1003) cx+cos(ang)*ri,cy+sin(ang)*ri,rdot,0,360.0
          write (iout,*) 'sk'
        end if
        riprev=ri
        angprev=ang
      end do
      write (iout,*) 'sk'
      do i=1,n
        cosav=cosav+res(1,i,incrid+idial)
        sinav=sinav+res(2,i,incrid+idial)
      end do
      cv=1.d0-dsqrt(cosav**2+sinav**2)/dfloat(n)
      lthick=2
      if (ndprow .gt. 6) lthick=1
      if (n .gt. 0) then
c       Draw average line
        call rgbcolor(iout,-1)
        write (iout,1004) lthick
        cosav=cosav/n
        sinav=sinav/n
        sqsum=sqrt(cosav**2+sinav**2)
        cosav=cosav/sqsum
        sinav=sinav/sqsum
        write (iout,*) 'np'
        write (iout,1000) cx+cosav*r0,cy+sinav*r0
        write (iout,1001) cx+cosav*(r0+r1),cy+sinav*(r0+r1)
        write (iout,*) 'sk'
c       Draw tick outside at the final position
        call rgbcolor(iout,-6)
        write (iout,1004) 2*lthick
        write (iout,*) 'np'
        rx=r*res(1,n,incrid+idial)
        ry=r*res(2,n,incrid+idial)
        write (iout,1000) cx+rx,cy+ry
        write (iout,1001) cx+1.1*rx,cy+1.1*ry
        write (iout,*) 'sk'
      end if
c     Draw CV bar
      call rgbcolor(iout,-9)
      write (iout,1004) 4
      write (iout,*) 'np'
      write (iout,1000) cx,cy+r0
      write (iout,1001) cx,cy+r0+r1*cv
      write (iout,*) 'sk'
      call rgbcolor(iout,9)
      write (iout,1004) 1
      if (ioutpr .gt. 0) then
        angmaxd=angmax*rdtodg
        angmind=angmin*rdtodg
        grd=ipdfgrd
        write (ioutpr,1006) lab(1:llab),angmind,angmaxd
        call zeroiti(npdf,0,360)
        call zeroit(pdf,360)
        do i=1,n
          ang=dacoscheck(ddd,res(1,i,incrid+idial),0,ioutpr,'DRAWDIAL')
          if (res(2,i,incrid+idial) .lt. 0.0) ang=pi2-ang
          ang=ang*rdtodg
c         ix=(ang-angmind)/grd+1
c         if (ix .gt. 360) ix=360
          ix=ang/grd+1
          npdf(ix)=npdf(ix)+1
        end do
c       ixmax=(angmaxd-angmind)/grd+1
c       if (ixmax .gt. 360) ixmax=360
        ixmin=angmind/grd+1
        ixmax=angmaxd/grd+1
        ixmax=min0(ixmax,360/ipdfgrd)
        pdfmax=0.0
c       write (88,*) 'R1=',r1,' ANGMIND,MAXD=',angmind,angmaxd
c       write (88,*) 'ANGMIN,MAX=',angmin,angmax
        do ix=ixmin,ixmax
          pdf(ix)=float(npdf(ix))/float(n)
          write (ioutpr,1007) ix,ix*grd-grd/2.0,pdf(ix)
          if (pdfmax .lt. pdf(ix)) pdfmax=pdf(ix)
        end do
        if (pdfmax .gt. 0.0 .and. mappdf .gt. 0) then
          do ix=ixmin,ixmax
            pdf(ix)=r0+r1+(pdf(ix)/pdfmax)*r1/2.75
          end do
          write (iout,1008) 'newpath'
          angmin2=angmin+(grd/2.0)/rdtodg
          write (iout,1000) cx+pdf(ixmin)*cos(angmin2),
     -      cy+pdf(ixmin)*sin(angmin2)
          do ix=ixmin+1,ixmax
            ang=(ix*grd-grd/2.0)/rdtodg
c           write (88,*) ix,' ANG=',ang
            write (iout,1001) cx+pdf(ix)*cos(ang),cy+pdf(ix)*sin(ang)
          end do
          if (ixmax .eq. 360/ipdfgrd .and. pdf(1) .gt. 0.0) then
            ang=(grd/2.0)/rdtodg
            write (iout,1001) cx+pdf(1)*cos(ang),cy+pdf(1)*sin(ang)
          end if
          write (iout,1008) 'stroke'
        end if
      end if
      return
1000  format(2f12.5,' m')
1001  format(2f12.5,' l')
1003  format(5f12.4,' arc')
1004  format(i5,' setlinewidth')
1005  format('/',a,' findfont',/,i2,' scalefont',/,'setfont')
1006  format(/,' Distribution of ',a,' in the range [',f7.2,',',f7.2,
     -  '] (deg)')
1007  format(i4,' ang=',f7.2,' P=',f7.5)
1008  format(a)
      end
