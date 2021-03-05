      subroutine rrconn(c,n,numres1,numres2,ifres,ilres,iatnum,
     -  ignoreh,irepuse,iadjtyp,iscalesum,resdistlim,nexpmax,npint,
     -  ipspage,npspages,ilog,iplot,line,index,ir1,ic1,
     -  imarkres,marks,imarks,iarep,iconnsum,itemp,rconnsum,inpfile,
     -  linpfile,markfile,lmarkfile,maxrec)
      dimension c(3,n),index(n),ifres(numres2),ilres(numres2),
     -  iatnum(n),imarks(numres2),iarep(n),iconnsum(n),itemp(n),
     -  rconnsum(n)
      character*1 marks(9)
      character* 132 line(maxrec)
      character*(*) inpfile,markfile
      parameter (MAXPHI=400,MAXCONN=2900,MAXCONN10=10*MAXCONN)
      parameter (IFILL8=MAXPHI*MAXPHI*MAXPHI-3*MAXCONN*MAXCONN
     -  -11*MAXCONN)
      common /nnwork/ rij1(MAXCONN,MAXCONN),rij2(MAXCONN,MAXCONN),
     -  rij3(MAXCONN,MAXCONN),xres(MAXCONN),yres(MAXCONN10),fill(IFILL8)
      dimension ifg(10),imf(10),iml(10),iexpplot(10),lfclab(10)
      character*8 fclab(10)
      character*200 title1,title2
      data connmax /0.0/
c     print *,'RRCONN numres1,numres2,nexpmax,iscalesum=',
c    -  numres1,numres2,nexpmax,iscalesum
c     print *,'RRCONN npspages,ipspage,maxrec=',
c    -  npspages,ipspage,maxrec
      write (6,1003) numres1,numres2
      if (iadjtyp .eq. 1)
     -   write (ilog,1005) '0/1 (integer)','matrix product'
      if (iadjtyp .eq. 2)
     -   write (ilog,1005) '0.0-1.0 (weighted)','matrix product'
      if (iadjtyp .eq. 3)
     -  write (ilog,1005) '0.0-1.0 (weighted)','matrix product-like sum'
      resdistlim2=resdistlim**2
      numres=numres2-numres1+1
      do i1=1,numres
        do i2=1,numres
          rij1(i1,i2)=0.0
          rij2(i1,i2)=0.0
          rij3(i1,i2)=0.0
        end do
      end do
      if (irepuse .gt. 0) then
        write (ilog,1001) 'representative atoms',resdistlim
        title2='Distances based on representative atoms;'
        ltitle2=40
c       Distances based on representative atoms
        do irr=numres1,numres2
c         First find the representative atom and atom range for residue irr
          call findat(iarep(irr),ifres(irr),ilres(irr),line,index,
     -     ir1,ic1,maxrec)
        end do
        do i1=numres1+1,numres2
          ii1=i1-numres1+1
          do i2=numres1,i1-1
            ii2=i2-numres1+1
            rijsq=dist2(c(1,iarep(i1)),c(1,iarep(i2)))
            if (rijsq .lt. resdistlim2) then
              rij1(ii1,ii2)=rijsq
              rij1(ii2,ii1)=rijsq
            end if
          end do
        end do
      else
c       Distances based on closest approach
        write (ilog,1001) 'closest approach',resdistlim
        title2='Distances based on closest approach;'
        ltitle2=36
        if (ignoreh .eq. 1) write (ilog,1002)
        do i1=numres1+1,numres2
          ii1=i1-numres1+1
          do i2=numres1,i1-1
            ii2=i2-numres1+1
            call findapproach(c,ifres(i1),ilres(i1),ifres(i2),
     -        ilres(i2),iatnum,ignoreh,iarep1,iarep2,rijsq,maxrec)
            if (rijsq .lt. resdistlim2) then
              rij1(ii1,ii2)=rijsq
              rij1(ii2,ii1)=rijsq
            end if
          end do
        end do
      end if
c     Convert distances to adjacencies
      do i1=2,numres
        do i2=1,i1-1
          if (rij1(i1,i2) .gt. 0.0) then
            if (iadjtyp .eq. 1) then
              rij=1.0
            else
              if (rij1(i1,i2) .lt. 9.0) then
                rij=1.0
              else
                r=sqrt(rij1(i1,i2))
                rij=(resdistlim-r)/(resdistlim-3.0)
              end if
            end if
            rij1(i1,i2)=rij
            rij1(i2,i1)=rij
            rij2(i1,i2)=rij
            rij2(i2,i1)=rij
            rij3(i1,i2)=rij
            rij3(i2,i1)=rij
          end if
        end do
      end do
      do i=1,numres
        do j=1,numres
          itemp(j)=rij1(i,j)
        end do
        write (ilog,1009) i,(itemp(j),j=1,numres)
      end do
c     Calculate the adjaceny matrix powers
      do i=1,10
        ifg(i)=1
        imf(i)=(i-1)*numres+1
        iml(i)=i*numres
      end do
      do i=1,numres
        xres(i)=i
      end do
      nexpplot=0
      do nexp=1,nexpmax
        if (nexp .gt. 1) then
          do i1=1,numres
            do i2=1,numres
              sum=0.0
              if (iadjtyp .lt. 3) then
                do k=1,numres
                  sum=sum+rij1(i1,k)*rij2(k,i2)
                end do
              else
                do k=1,numres
                  if (rij1(i1,k) .gt. 0.0 .and. rij2(k,i2) .gt. 0.0)
     -              sum=sum+rij1(i1,k)+rij2(k,i2)
                end do
              end if
              rij3(i1,i2)=sum
            end do
          end do
        end if
        do i1=1,numres
          colsum=0.0
          do i2=1,numres
            colsum=colsum+rij3(i1,i2)
            rij2(i1,i2)=rij3(i1,i2)
          end do
          rconnsum(i1)=colsum
        end do
        connmax=0.0
        do i=1,numres
          if (rconnsum(i) .gt. connmax) connmax=rconnsum(i)
        end do
        if (iscalesum .eq. 1) then
          do i=1,numres
            iconnsum(i)=100.0*rconnsum(i)/connmax
          end do
        else
          scalefac=1.0
          nlog=alog10(connmax)
          if (nlog .gt. 3) scalefac=10.0**(2-nlog)
          do i=1,numres
            iconnsum(i)=rconnsum(i)*scalefac
          end do
        end if
        write (ilog,1000) nexp,(iconnsum(i),i=1,numres)
        if ((mod(nexp,npint) .eq. 1 .or. nexp .eq. nexpmax .or.
     -      npint .eq. 1) .and. nexpplot .lt. 10) then
          nexpplot=nexpplot+1
          iexpplot(nexpplot)=nexp
          write (fclab(nexpplot),1004) nexp
          lfclab(nexpplot)=8
          if (iscalesum .eq. 0) then
            do i=1,numres
              yres(imf(nexpplot)-1+i)=iconnsum(i)
            end do
            ny=5
            ydiv=0.0
          else
            do i=1,numres
              yres(imf(nexpplot)-1+i)=(nexpplot-1)+rconnsum(i)/connmax
            end do
            ny=nexpplot
            ydiv=1.0
          end if
        end if
      end do
c     Calculate correlation between marks and high/low col sums
      nmarksum=0
      nmarkup=0
      nmarkdown=0
      nsumup=0
      nsumdown=0
      do i=1,numres
        if (rconnsum(i) .gt. connmax/2.0) then
          nsumup=nsumup+1
        else
          nsumdown=nsumdown+1
        end if
        if (imarks(i) .gt. 0) then
          nmarksum=nmarksum+1
          if (rconnsum(i) .gt. connmax/2.0) then
            nmarkup=nmarkup+1
          else
            nmarkdown=nmarkdown+1
          end if
        end if
      end do
      if (nmarksum .gt. 0) then
        write (ilog,1007) nmarksum,float(nmarkup)/float(nmarksum),
     -    float(nmarkdown)/float(nmarksum),numres,
     -    float(nsumup)/float(numres),float(nsumdown)/float(numres)
      else
        write (ilog,1008) inpfile(1:linpfile)
      end if
      title1='Adjacency matrix analysis for file '
      ltitle1=35
      title1(ltitle1+1:ltitle1+linpfile)=inpfile(1:linpfile)
      ltitle1=ltitle1+linpfile
      write (title2(ltitle2+1:ltitle2+14),1006) resdistlim
      ltitle2=ltitle2+14
      if (iadjtyp .eq. 1) then
        title2(ltitle2+1:ltitle2+27)=' 0/1 matrix, matrix product'
        ltitle2=ltitle2+27
      else
        title2(ltitle2+1:ltitle2+34)=' 0.0-1.0 matrix, product-like sum'
        ltitle2=ltitle2+34
      end if
      if (imarkres .gt. 0) then
        title2(ltitle2+1:ltitle2+13)='; Mark file: '
        ltitle2=ltitle2+13
        title2(ltitle2+1:ltitle2+lmarkfile)=markfile(1:lmarkfile)
        ltitle2=ltitle2+lmarkfile
      end if
      call rounddiv(numres,10,nx,nxdiv)
      xdiv=nxdiv
      call plotnps(xres,yres,MAXCONN,MAXCONN10,nexpplot,imf,iml,ifg,
     -  0.0,1.0,0.0,xdiv,nx,0.0,ydiv,ny,1,ltitle1,title1,ltitle2,title2,
     -  'Residue #',9,fclab,lfclab,imarkres,imarks,marks,iplot,ipspage,
     -  npspages,inperr,ilog)
      return
1000  format(' Bond density (scaled); Adjacency matrix power=',i2,/,
     -  (20i4))
1001  format(' Adjacency is based on ',a,/,
     -  ' Distance threshold =',f5.2,' A')
1002  format(' Closest approach is based on heavy atoms only')
1003  format(' Adjacency matrix analysis for the residue range [',i5,
     -  ',',i5,']')
1004  format('Power=',i2)
1005  format(' Matrix: ',a,' Operation: ',a)
1006  format(' RRmax=',f4.1,' A;')
1007  format(' # of marks=',i4,' x up=',f5.3,' x down=',f5.3,
     -  ' # of res=',i4,' x up=',f5.3,' x down=',f5.3)
1008  format(' No marks were found for file ',a)
1009  format(i4,1x,10i1,1x,10i1,1x,10i1,1x,10i1,1x,10i1,1x,/,
     -  (5x,5(10i1,1x)))
      end
