      subroutine resautocorr(ires,incrementac,ncf,xyplot,mxframes)
      dimension xyplot(2,mxframes)
      parameter (MAXFRAMES=50000,MAXCOPY=600)
      common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
     -  increment,increment2,res(2,MAXFRAMES,MAXCOPY),
     -  x0(MAXCOPY),y0(MAXCOPY),nxselres,ixselres(MAXCOPY)
c     Calculate the relaxation time of the ires-th stored residue
c     print *,'RESAUTOCORR ires,nframe,increment=',ires,nframe,increment
      ncf=nframe/(2*incrementac)
      if (ncf .eq. 0) then
        write (6,1000) nframe,incrementac
        return
      end if
      call zeroit(xyplot,2*ncf)
      do i=1,ncf
        i0=i*incrementac
        do inc=1,ncf-i
          i1=(i+inc)*incrementac
          xyplot(1,inc+1)=xyplot(1,inc+1)+
     -      (res(1,i0,2*ires-1)*res(1,i1,2*ires-1)+
     -      res(2,i0,2*ires-1)*res(2,i1,2*ires-1)+
     -      res(1,i0,2*ires)*res(1,i1,2*ires)+
     -      res(2,i0,2*ires)*res(2,i1,2*ires))/2.0
        end do
      end do
      xyplot(1,1)=1.0
      do i=2,ncf-1
        xyplot(1,i)=xyplot(1,i)/float(ncf-i)
      end do
      xyplot(1,ncf)=xyplot(1,ncf-1)
      return
1000  format(' PROGRAM ERROR in resautocorr: nframe (',i8,
     -  ') < 2* increment (',i8,')')
      end
