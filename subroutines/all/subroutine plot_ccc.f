      subroutine plot_ccc(ndials,iout,ips,ipspage,angnames,langnames,
     -  inpfile,namleni)
      dimension langnames(ndials)
      character*(*) angnames(ndials)
      character*(*) inpfile
      common /colorinfo/ ncolcode,maxcolcode
      parameter (MAXPHI=400,MAX2D=5000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-(MAX2D*MAX2D+MAX2D))
      common /nnwork/ ccc(MAX2D,MAX2D),itemp(MAX2D),fill(IFILL2)
      dimension kc(1,1),xyval(1)
      real*8 dc(1,1)
      character*4 yclab(1)
      character*80 title2
      data nyclab /1/,lyclab /1/
      data title2 /'Circular correlation matrix'/,ltitle2 /28/
      if (ips .eq. 0) return
      rcmin=-1.0
      rcmax=1.0
      call getint('Number of colors to use',23,5,1,maxcolcode,
     -  ncolcode,0)
      scalefac=amin1(1.0,500.0/float(ndials))
      ixdel=25
!     iydel=115
      iydel=100
      iytop=0
      incinp=max0(1,500/ndials)
      call indexit(itemp,1,ndials,0)
      call plotmat(ips,kc,ccc,dc,ndials,ndials,0,0,0,0,1,0,0,iydel,00,
     -  iytop,rcmin,rcmax,ncolcode,maxcolcode,ixdel,iydel,incinp,
     -  scalefac,itemp,itemp,itemp,inpfile,namleni,title2,ltitle2,0,' ',
     -  0,xyval,yclab,nyclab,lyclab,1,MAX2D,1,MAX2D,MAX2D,ipspage,0)
      ixd=ixdel
      if (ncolcode .le. 5) ixd=ixd+65
      call colcodeminmax(ips,ixd,-65,0,ncolcode,maxcolcode,rcmin,
     -  rcmax)
      call pswrite(ips,150,25,'m',1)
      write (ips,*) '(Range of data: [-1.0,1.0]) show'
      return
      end
