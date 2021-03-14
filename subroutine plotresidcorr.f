      subroutine plotresidcorr(ncorr,nframe,index,indexr,
     -  ncol,maxcol,nrep,iout,iplt,xm,ym,title,ltitle,
     -  ipspage,iucorrmat,icovmatplot,temp,maxtemp)
      dimension index(ncorr),indexr(ncorr)
      character*(*) title
      dimension temp(maxtemp)
      parameter (MAXPHI=400,MAX2D=5000)
      parameter (IFILL4=MAXPHI*MAXPHI*MAXPHI-
     -  (2*MAX2D*MAX2D+17*MAX2D))
      real*8 trajcorr,cav1,cav2,cavs1,cavs2
      common /nnwork/ trajcorr(MAX2D,MAX2D),cav1(3,MAX2D),
     -  cav2(3,MAX2D),cavs1(MAX2D),cavs2(MAX2D),
     -  row(MAX2D),fill(IFILL4)
      real*8 ddot
      dimension mx(1,1),rx(1,1),ixshuffle(MAX2D)
      character*4 yclab(1)
      data nyclab /1/,lyclab /1/
      call indexit(ixshuffle,1,MAX2D,0)
      call write_traj_lim(iout,'Residue covariances and correlations',
     -  36,1,incr_tr,0)
      write (iout,1000) ncorr,nframe,(index(ir),ir=1,ncorr)
      do ir=1,ncorr
        do k=1,3
          cav1(k,ir)=cav1(k,ir)/nframe
          cav2(k,ir)=cav2(k,ir)/nframe
        end do
        cavs1(ir)=cavs1(ir)/nframe
        cavs2(ir)=cavs2(ir)/nframe
      end do
      rcmax=-100000.0
      rcmin=-rcmax
      do ir=1,ncorr
        do jr=1,ncorr
          trajcorr(ir,jr)=(trajcorr(ir,jr)/nframe-
     -      ddot(cav1(1,ir),cav2(1,jr)))
          if (trajcorr(ir,jr) .lt. rcmin) rcmin=trajcorr(ir,jr)
          if (trajcorr(ir,jr) .gt. rcmax) rcmax=trajcorr(ir,jr)
        end do
      end do
      write (iout,1004) rcmin,rcmax
      do ir=1,ncorr
        if (rcmin .lt. 10000.0 .or. rcmax .gt. 10000.0) then
          write (iout,1003) (trajcorr(ir,jr),jr=1,ncorr)
        else
          write (iout,1002) (trajcorr(ir,jr),jr=1,ncorr)
        end if
      end do
c     Save the covariance matrix for eigenvalue calculation
      if (iucorrmat .gt. 0) then
        call openfile(iucorrmat,0,'log',3,'new','trajcorr.mat',12,
     -      notfnd,0,2,1,0,0)
        do jr=1,ncorr
          write (iucorrmat) (trajcorr(ir,jr),ir=1,ncorr)
        end do
      end if
      inc=max0(1,500/ncorr)
      scalefac=amin1(1.0,500.0/float(ncorr))
      iydel=150
      iytop=0
      if (icovmatplot .eq. 1) then
        call plotmat(iplt,mx,rx,trajcorr,ncorr,ncorr,0,0,0,0,1,nrep,30,
     -    iydel,0,iytop,rcmin,rcmax,ncol,maxcol,ixdelsh,iydelsh,inc,
     -    scalefac,indexr,ixshuffle,ixshuffle,title,ltitle,' ',0,1,' ',
     -    0,temp,yclab,nyclab,lyclab,1,1,MAX2D,MAX2D,MAX2D,
     -    ipspage,0)
        iydel=iydel-50
        call rainbowscale(iplt,ixdelsh+50,450,iydel,0,0.0,rcmin,rcmax,
     -    'Covariance',10)
        call plothead(iplt,xm,ym-15,title,0,
     -    'Residue covariance matrix',25,' ',0,' ',0)
        write (iplt,*) 'showpage'
      end if
      do ir=1,ncorr
        cavs1(ir)=dsqrt(dabs(trajcorr(ir,ir)))
      end do
      do ir=1,ncorr
        do jr=1,ir-1
          trajcorr(ir,jr)=trajcorr(ir,jr)/(cavs1(ir)*cavs1(jr))
          trajcorr(jr,ir)=trajcorr(jr,ir)/(cavs1(ir)*cavs1(jr))
        end do
      end do
      do ir=1,ncorr
        trajcorr(ir,ir)=1.d0
      end do
      write (iout,*) 'The correlation matrix:'
      do ir=1,ncorr
        write (iout,1001) (trajcorr(ir,jr),jr=1,ncorr)
      end do
      iydel=150
      iytop=0
      call plotmat(iplt,mx,rx,trajcorr,ncorr,ncorr,0,0,0,0,1,nrep,30,
     -  iydel,0,iytop,-1.0,1.0,ncol,maxcol,ixdelsh,iydelsh,inc,scalefac,
     -  indexr,ixshuffle,ixshuffle,title,ltitle,' ',0,1,' ',0,temp,
     -  yclab,nyclab,lyclab,1,1,MAX2D,MAX2D,MAX2D,ipspage,0)
      iydel=iydel-50
      call rainbowscale(iplt,ixdelsh+50,450,iydel,0,0.0,-1.0,1.0,
     -  'Correlation',11)
      call plothead(iplt,xm,ym-15,title,0,
     -  'Residue correlation matrix',26,' ',0,' ',0)
      write (iplt,*) 'showpage'
      return
1000  format(' Number of residues used=',i5,/,
     -  ' Number of structures used=',i6,/,
     -  ' Atomindices used for the calculation:',/,(10i6))
1001  format(10f8.4)
1002  format(10f8.2)
1003  format(5e12.5)
1004  format(' The covariance matrix (min=',e12.5,' max=',e12.5,'):')
      end
