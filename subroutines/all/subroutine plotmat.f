      subroutine plotmat(ips,kc,rc,dc,nresx,nresy,ixinc,iyinc,
     -  ixincshift,iyincshift,navg,nrep,ixdel0,iydel0,idown,iytop,rcmin,
     -  rcmax,ncol,maxcol,ixdel,iydel,incinp,scalefacinp,index,
     -  ixshufflex,ixshuffley,title,ltitle,title2,ltitle2,itrajname,
     -  xylab,lxylab,xyval,yclab,nyclab,lyclab,mxr,mxrr,mxrd,mx,mxix,
     -  ipspage,ipsclose)
      character*(*) title,title2,xylab,yclab(nyclab)
      real*8 dc
      dimension kc(mxr,mxr),rc(mxrr,mxrr),dc(mxrd,mxrd),index(mxix),
     -  ixshufflex(mx),ixshuffley(mx),xyval(mx)
      data ix0 /0/,iy0 /0/,rcxy /0.0/
c     kc,rc,dc: matrix in integer, real*4 and real*8 formats
c     mxr,mxrr,mxrd: dimensions of kc,rc,dc (only one of them can be > 1)
c     nresx, nresy: matrix dimension
c     ishufflex, ixinc: X axis label of the ith row is ixshufflex(i+ixinc)
c     ishuffley, iyinc: Y axis label of the ith row is ixshuffley(i+iyinc)
c     ixincshift,iyincshift: deduct from the matrix indices
c     navg: number of residues the matrix entries were averaged over
c     nrep: When >1 then it was a repeat plot
c     ixdel,iydel: increments in the x and y direction
c     scalefac: Overall scaling factor (500*550 matrix fills the page
c     with scalefac=1
c     rcmin,rcmax: min and max of matrix value allowed
c     ncol, maxcol: requested and maximum number of colors
c     inc: increment per matrix row/col
c     print *,'PLOTMAT ips,ncol=',ips,ncol,' nresx,nresy=',nresx,nresy
c     print *,'PLOTMAT mxr,mxrr,mxrd=',mxr,mxrr,mxrd
c     print *,'PLOTMAT mx,mxix=',mx,mxix
c     print *,'PLOTMAT iydel,iydel0,idown,iyinc=',
c    -  iydel,ixdel0,idown,iyinc
c     print *,'PLOTMAT ixdel,ixdel0,iytop=',iytop,ixdel,ixdel0
c     print *,'PLOTMAT ltitle,ltitle2=',ltitle,ltitle2,
c    -  ' xylab=',xylab(1:lxylab)
c      do i=1,10
c        write (6,7733) i,(kc(i,j),j=1,10)
c      end do
c7733  format(' KC i=',i3,' kc(i)=',10i3)
c     if (ltitle .gt. 0) print *,'PLOTMAT title=',title(1:ltitle)
c     if (ltitle2 .gt. 0) print *,'PLOTMAT title=',title2(1:ltitle2)
c     print *,'PLOTMAT ITRAJNAME=',itrajname
      incr_tr=1
      incr_tr2=1
      scalefac=scalefacinp
      lycinc=lyclab/10
      scalefac=scalefac*(1.0-lycinc*0.04)
      icharsize=12
      ipsdraw=ips
      ipsdraw=-ipsdraw
      maxrrgb=maxcol
      maxcol1=maxcol+1
c     write (6,8422) 'X',(ixshufflex(i),i=1,nresx)
c     write (6,8422) 'Y',(ixshuffley(i),i=1,nresy)
c8422 format(' IXSHUFFLE',a1,':',/,(20i4))
      ipspage=ipspage+1
      write (ips,1020) ipspage
      inc=incinp*navg
      write (ips,3011) icharsize
      if (rcmax .eq. rcmin) rcmax=rcmax+0.1
      ixdel=ixdel0
      ixdelmat=ixdel0
      iydel=iydel0
      if (scalefac .ne. 1.0) then
        write (ips,3003) scalefac,scalefac
        icharscale=float(icharsize)/scalefac
        write (ips,3011) icharscale
        iydel=iydel/scalefac
        iytop=iytop/scalefac
        ixdel=ixdel/scalefac
      end if
      if (nresx*scalefac .lt. 550) ixdelmat=ixdel0+40/scalefac
      if (iytop .eq. 0) iytop=iydel+nresy*inc
      if (ltitle .gt. 0) then
        iytop=iytop+15
        call pswrite(ips,ixdel,iytop,'m',1)
        write (ips,*) '(System: ) show'
        call psshow(ips,title,ltitle)
      end if
      if (itrajname .eq. 3) then
c       Plot both names
        iytop=iytop+15
        call pswrite(ips,ixdel,iytop,'m',1)
        call write_traj_lim(ips,' ',1,2,incr_tr2,1)
        iytop=iytop+15
        call pswrite(ips,ixdel,iytop,'m',1)
        call write_traj_lim(ips,' ',1,1,incr_tr,1)
      else if (itrajname .gt. 0) then
        iytop=iytop+15
        call pswrite(ips,ixdel,iytop,'m',1)
        call write_traj_lim(ips,' ',1,itrajname,incr_tr,1)
      end if
      if (ltitle2 .gt. 0) then
        iytop=iytop+15
        call pswrite(ips,ixdel,iytop,'m',1)
        call psshow(ips,title2,ltitle2)
      end if
c     print *,'PLOTMAT iydel,ixinc,iyinc,inc=',
c    -  iydel,ixinc,iyinc,inc
      write (ips,3005) 'Plot matrix elements'
      do iyy=1,nresy
        iy=ixshuffley(iyy+iyinc)-iyincshift
c       write (77,*) 'iyy,iyinc,iy=',iyy,iyinc,iy
        kcprev=-1
        do ixx=1,nresx
          ix=ixshufflex(ixx+ixinc)-ixincshift
          if (mxr .gt. 1) then
            kccurr=kc(ix,iy)
          else
c           Convert rc(ix,iy) or dc(ix,iy) into color code
            if (mxrr .gt. 1) rcxy=rc(ix,iy)
            if (mxrd .gt. 1) rcxy=dc(ix,iy)
            if (ncol .eq. 0) then
              kccurr=100.0*(rcxy-rcmin)/(rcmax-rcmin)+1.0
            else
              if (rcxy .le. rcmin) then
                kccurr=0
              else if (rcxy .gt. rcmax) then
                kccurr=ncol+1
              else
                rr=amax1(-0.00001,rcxy-rcmin-0.00001)
                kccurr=(rr/(rcmax-rcmin))*ncol+1
              end if
            end if
          end if
c          write (77,7688) ix,iy,rc(ix,iy),kccurr
c7688      format('ix,iy=',2i6,' rc=',f10.5,' kcurr=',i3)
          if (nrep .le. 1) then
            if (kccurr .ne. kcprev) then
              if (kcprev .ne. -1 .and. (kcprev .ne. maxcol1 .or.
     -            ncol .eq. 0)) then
c               Close the rectangle that is open
                call pswrite(ips,ixdelmat+(ixx-1)*inc,iydel+(iyy)*inc,
     -            'l',1)
                call pswrite(ips,ixdelmat+(ixx-1)*inc,iydel+(iyy-1)*inc,
     -            'l',1)
                call pswrite(ips,ix0-inc,iy0-inc,'l',1)
                write (ips,3006)
                ix0=-1
              end if
              if (kccurr .ne. maxcol1 .or. ncol .eq. 0) then
c               Open new rectangle
                if (ncol .eq. 0) then
                  call rrgbcolor(ips,kccurr,100,0)
                else
                  call rgbcolor(ips,kccurr)
                end if
                ix0=ixdelmat+(ixx)*inc
                iy0=iydel+(iyy)*inc
c                write (77,8766) ixx,iyy,ix0,iy0,inc
c8766  format(' ixx,iyy=',2i4,' ix0,iy0=',2i4,' inc=',i2)
                call pswrite(ips,ix0-inc,iy0-inc,'m',1)
                call pswrite(ips,ix0-inc,iy0,'l',1)
              end if
              kcprev=kccurr
            end if
          end if
        end do
        if (ix0 .ne. -1 .and. nrep .le. 1) then
c         Close last opened rectangle
          call pswrite(ips,ixdelmat+(nresx)*inc,iydel+(iyy)*inc,'l',1)
          call pswrite(ips,ixdelmat+(nresx)*inc,iydel+(iyy-1)*inc,'l',1)
          call pswrite(ips,ix0-inc,iy0-inc,'l',1)
          write (ips,3006)
        end if
      end do
c     Write residue/frame number scale
c     print *,'Matrix entries plotted'
c     Plot protein id
      ridown=idown
      idown1=idown+06
      idown2=idown1+9
      rydel=iydel-idown
      rydel1=iydel-idown1
      rydel2=iydel-idown2
c     if (nrep .le. 1) then
c       call rgbcolor(ips,9)
c       call pswrite(ips,ixdel,iydel-idown1,'m',1)
c     end if
      write (ips,3004) 2
      call rgbcolor(ips,9)
      irxmax=nresx
      irymax=nresy
      if (lxylab .gt. 1) then
        write (ips,3000)
        ix=ixdelmat+(nresx*inc)*0.45
        iydown=30/scalefac
        call pswrite(ips,ix,iydel-idown-iydown,'m',1)
        write (ips,3009) xylab(1:lxylab)
        write (ips,3001)
        xmax=xyval(nresx)
        call roundlim(xmax,xdiv,nxdiv)
        ymax=xyval(nresy)
        if (iyinc .gt. 0) ymax=xyval(nresy+iyinc)-xyval(iyinc)
        call roundlim(ymax,ydiv,nydiv)
        xdiv_i=xdiv*float(nresx)/xmax
        ydiv_i=ydiv*float(nresy)/ymax
        iydiv=ydiv
      else
        call roundlimint(nresx,ixdiv,nxdiv)
        call roundlimint(nresy,iydiv,nydiv)
        xdiv=ixdiv
        ydiv=iydiv
        xdiv_i=ixdiv
        ydiv_i=iydiv
        xmax=xdiv*nresx
        ymax=ydiv*nresy
      end if
c     print *,'DRAWRECT nresy,iydiv,iydel=',nresy,iydiv,iydel
c     print *,'DRAWRECT ixdelmat,inc=',ixdelmat,inc
      write (ips,3005) 'Plot bounding rectangle'
      call drawrect(ipsdraw,9,1,ixdelmat,ixdelmat+nresx*inc,iydel,
     -  iydel+max0(nresy,iydiv/incr_tr)*inc,nrep)
      write (ips,3005) 'Plot ticks, axis values'
      ixadd=0
      iyadd=0
      if (abs(xmax-nxdiv*xdiv) .lt. xmax*0.01) ixadd=1
      if (abs(ymax-nydiv*ydiv) .lt. ymax*0.01) iyadd=1
      idown1=idown1/scalefac
      idown2=idown2/scalefac
      ixsshift=14/scalefac
      itickstep=10/(xdiv_i*inc)+1
      do itic=1,nxdiv+ixadd,itickstep
        if (lxylab .le. 1) then
          ir=itic*xdiv_i*inc
        else
          ir=(itic-1)*xdiv_i*inc
        end if
        irinc=max0(1,ir/inc)
c       write (6,8768) nresx,irinc,xdiv,xdiv_i,index(nresx)
c8768   format(' nresx=',i6,' irinc=',i6,' xdiv,xdiv_i=',2f6.2,
c    -    ' index(nresx)=',i6)
        if (index(irinc) .le. irxmax) then
          if (nrep .le. 1) then
c           Upper tick
            write (ips,3000)
            call pswrite(ips,ixdelmat+ir,iydel+nresy*inc+idown,'m',1)
            call pswrite(ips,ixdelmat+ir,iydel+nresy*inc+idown1,'l',1)
            write (ips,3001)
c           Lower tick
            write (ips,3000)
            call pswrite(ips,ixdelmat+ir,iydel-idown,'m',1)
            call pswrite(ips,ixdelmat+ir,iydel-idown1,'l',1)
            write (ips,3001)
c           Scale
            call pswrite(ips,ixdelmat+ir-ixsshift,iydel-idown2,'m',1)
            if (lxylab .le. 1) then
              write (ips,3008) index(irinc)
            else if (xylab(1:1) .eq. ' ' .or.
     -               xdiv .ge. 1.0) then
c    -               xdiv .ge. 100.0) then
              ix=navg*(itic-1)*xdiv
              write (ips,3008) ix
            else
              write (ips,3010) navg*(itic-1)*xdiv
            end if
          end if
        end if
      end do
      itickstep=10/(ydiv_i*inc)+1
      do itic=1,nydiv+iyadd,itickstep
        if (lxylab .le. 1) then
          ir=itic*ydiv_i*inc
        else
          ir=(itic-1)*ydiv_i*inc
        end if
        irinc=max0(1,ir/inc)+iyinc
        if (index(irinc) .le. iyinc+irymax .or. itic .eq. 1) then
c         Left tick
          write (ips,3000)
          call pswrite(ips,ixdelmat,iydel+ir,'m',1)
          call pswrite(ips,ixdelmat-(idown1-idown),iydel+ir,'l',1)
          write (ips,3001)
c         Right tick
          write (ips,3000)
          call pswrite(ips,ixdelmat+nresx*inc,iydel+ir,'m',1)
          if (nyclab .le. 1) then
            call pswrite(ips,ixdelmat+nresx*inc+(idown1-idown),iydel+ir,
     -        'l',1)
          else
            call pswrite(ips,ixdelmat+nresx*inc-(idown1-idown),iydel+ir,
     -        'l',1)
          end if
          write (ips,3001)
          if (ixdelmat .gt. ixdel0) then
c           Scale
            call pswrite(ips,ixdel0,iydel+ir-4,'m',1)
            if (lxylab .le. 1) then
              write (ips,3008) index(irinc)
            else if (xylab(1:1) .eq. ' ' .or.
     -               ydiv .ge. 1.0) then
c    -               ydiv .ge. 100.0) then
              iy=navg*(itic-1)*ydiv+iyinc
              write (ips,3008) iy
            else
              write (ips,3010) navg*(itic-1)*ydiv
            end if
          end if
        end if
      end do
c     print *,'PLOTMAT NYCLB,LYCLAB=',nyclab,lyclab
      if (nyclab .gt. 1) then
c       Print residue identifiers
        write (ips,3011) icharsize/2+1
        do iy=1,nyclab
          call pswrite(ips,ixdelmat+nresx*inc,iydel+(iy-1)*inc+inc/3,
     -      'm',1)
          write (ips,3009) yclab(iy)(1:lyclab)
        end do
        write (ips,3011) icharsize
      end if
      if (scalefac .ne. 1.0) then
        write (ips,3003) 1.0/scalefac,1.0/scalefac
        write (ips,3011) icharsize
      end if
      if (ipsclose .eq. 1) close (ips)
      return
1020  format('%%Page: 1 ',i4)
3000  format('np')
3001  format('sk')
3003  format(f10.6,f10.6,' scale')
3004  format(i5,' lw')
3005  format('% ',a)
3006  format('f')
3008  format('(',i6,') show')
3009  format('(',a,') show')
3010  format('(',f6.2,') show')
3011  format('/Helvetica findfont',/,i4,' scalefont',/,'setfont')
      end
