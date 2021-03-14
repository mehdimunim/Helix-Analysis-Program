      subroutine ionreset(nmolfst,nmolslt,nmolshift,c,molsltlim,it1,atw,
     -  cell,ncell,ixyzexcld,ixyzincld,maxrsd,mxat)
      dimension c(3,mxat),molsltlim(3,maxrsd),it1(mxat),atw(mxat),
     -  cell(3,ncell)
      dimension c0(3),cmin(3),cmax(3)
c     print *,'MOLRESET ncell,ixyzexcld,ixyzincld=',
c    -  ncell,ixyzexcld,ixyzincld
      nmolshift=0
      do is=nmolfst,nmolslt
        if (molsltlim(3,is) .eq. 0) then
          call extension(c,it1,0,molsltlim(1,is),molsltlim(2,is),
     -      cmin,cmax,c0,0,0,v)
        else if (molsltlim(3,is) .eq. -1) then
          call cofms(c(1,molsltlim(1,is)),c0,
     -      molsltlim(2,is)-molsltlim(1,is)+1,atw)
        else
          call trnsfr(c0,c(1,molsltlim(3,is)),3)
        end if
        call genimdist123dim(c0,cell,1,ncell,ixyzexcld,ixyzincld,
     -    img,rmin2)
        nrep=0
        do while (img .gt. 1 .and. nrep .lt. 10)
c         if (is .lt. 4) print *,'MOLRESET is=',is,' img=',img
          call shiftmol(c(1,molsltlim(1,is)),
     -      molsltlim(2,is)-molsltlim(1,is)+1,cell(1,img),
     -      c(1,molsltlim(1,is)),-1.0)
          if (img .gt. 1) nmolshift=nmolshift+1
c         See if image is now in the central cell
          call shiftmol(c0,1,cell(1,img),c0,-1.0)
          call genimdist123dim(c0,cell,1,ncell,ixyzexcld,ixyzincld,
     -      img,rmin2)
          nrep=nrep+1
        end do
      end do
      return
      end
