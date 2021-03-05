      subroutine plot_atomdist_sd(n,line,index,inamcol1,inamcol2,
     -  irescol1,irescol2,iresncol1,iresncol2,ldist,ndist,nframe,indexa,
     -  ixshuffle,xtraj,title,trajfile,ltrajfile,isdtyp,iw0,iw1,
     -  ipspage,maxrec)
      dimension ldist(ndist),index(n),indexa(ndist),
     -  ixshuffle(ndist),xtraj(ndist)
      character*132 line(maxrec)
      character*80 title
      character*(*) trajfile
      parameter (MAXPHI=400,MAX2D=5000)
      parameter (IFILL9=MAXPHI*MAXPHI*MAXPHI-
     -  (2*MAX2D*MAX2D+11*MAX2D))
      real*8 trajdist,cav,cavs,rr,dscprod
      common /nnwork/ trajdist(MAX2D,MAX2D),cav(3,MAX2D),
     -  cavsng(3,MAX2D),cavs(MAX2D),fill(IFILL9)
      common /colorinfo/ ncolcode,maxcolcode
      dimension kc(1,1),rc(1,1),lnormtyp(4)
      character*4 yclab(1)
      character*17 normtyp(4)
      data normtyp /' ','/RMSD(i)*RMSD(j) ','/r(i,j) ','*r(i,j) '/,
     -  lnormtyp /1,17,8,8/
c     print *,'PLOT ATOMDIST_SD nframe,ndist=',nframe,ndist
      data nyclab /1/,lyclab /1/
      rcmax=0.0
      rcmin=100000.0
      rmsdmax=0.0
      rmsdmin=100000.0
      do ia=1,ndist
        xtraj(ia)=ia
        do k=1,3
          cav(k,ia)=cav(k,ia)/dfloat(nframe)
          cavsng(k,ia)=cav(k,ia)
        end do
        rr=dscprod(cav(1,ia),cav(1,ia))
        cavs(ia)=dsqrt(cavs(ia)/dfloat(nframe)-rr)
        if (rmsdmax .lt. cavs(ia)) rmsdmax=cavs(ia)
        if (rmsdmin .gt. cavs(ia)) rmsdmin=cavs(ia)
        do ja=1,ia-1
          sd=dsqrt(trajdist(ia,ja)/dfloat(nframe)-
     -      (trajdist(ja,ia)/dfloat(nframe))**2)
          if (isdtyp .eq. 2) then
            sd=sd/(cavs(ia)*cavs(ja))
          else if (isdtyp .eq. 3) then
            sd=sd/sqrt(dist2(cavsng(1,ia),cavsng(1,ja)))
          else if (isdtyp .eq. 3) then
            sd=sd*sqrt(dist2(cavsng(1,ia),cavsng(1,ja)))
          end if
          trajdist(ia,ja)=sd
          trajdist(ja,ia)=sd
          if (rcmin .gt. sd) rcmin=sd
          if (rcmax .lt. sd) rcmax=sd
        end do
        trajdist(ia,ia)=0.d0
      end do
      do ia=1,ndist
        write (iw0,1000) ia,ldist(ia),
     -    line(index(ldist(ia)))(inamcol1:inamcol2),
     -    line(index(ldist(ia)))(irescol1:irescol2),
     -    line(index(ldist(ia)))(iresncol1:iresncol2),
     -    (cav(k,ia),k=1,3),cavs(ia),
     -    normtyp(isdtyp)(1:lnormtyp(isdtyp)),
     -    (trajdist(ia,ja),ja=1,ndist)
      end do
      write (6,1011) 'SD(rij)',normtyp(isdtyp)(1:lnormtyp(isdtyp)),
     -  rcmin,rcmax
      write (iw0,1011) 'SD(rij)',normtyp(isdtyp)(1:lnormtyp(isdtyp)),
     -  rcmin,rcmax
      write (6,1011) 'site RMSD',' ',rmsdmin,rmsdmax
      write (iw0,1011) 'site RMSD',' ',rmsdmin,rmsdmax
      if (rmsdmax .gt. 10.0) then
        write (6,1006)
        write (iw0,1006)
      end if
      if (rcmin .eq. 0.0) then
        write (6,1017) 'minimum','some'
        write (iw0,1017) 'minimum','some'
      end if
      if (rcmax .eq. 0.0) then
        write (6,1017) 'maximum','all'
        write (iw0,1017) 'maximum','all'
        print *,'Exiting plotting'
        return
      end if
      navg=1
      nrep=0
      ixdelsh=0
      inc=max0(1,500/ndist)
      scalefac=amin1(1.0,500.0/float(ndist))
      iyrange=scalefac*(ndist*inc)+20
      ym_2d=800.0
      iytop=ym_2d*0.83
      iydel=max0(150,iytop-iyrange)
      iymax=iydel+iyrange
      write (iw1,1005) 25,iymax
      write (iw1,1003) normtyp(isdtyp)(1:lnormtyp(isdtyp)),
     -  trajfile(1:ltrajfile)
      iymax=iymax+15
      write (iw1,1005) 25,iymax
      call plotmat(iw1,kc,rc,trajdist,ndist,ndist,0,0,0,0,navg,
     -  nrep,25,iydel,00,iytop,rcmin,rcmax,ncolcode,maxcolcode,ixdelsh,
     -  iydelsh,inc,scalefac,indexa,ixshuffle,ixshuffle,title,0,
     -  ' ',0,0,'Atom #',6,xtraj,yclab,nyclab,lyclab,MAX2D,
     -    MAX2D,MAX2D,1,1,ipspage,0)
      if (ndist .ge. 50) iydel=iydel-50
      if (ndist .lt. 50) iydel=iydel-40
      ixcent=amax1(0.0,(scalefac*(ndist*inc)-80*ncolcode)/2)
      call colcodeminmax(iw1,20+ixcent,-iydel,nrep,ncolcode,
     -  maxcolcode,rcmin,rcmax)
      iydel=iydel-40
      write (iw1,1005) 50,iydel
      write (iw1,1002) rcmin,rcmax
      write (iw1,1005) 50+250,iydel
      return
1000  format(i5,' Atom',i7,' (',a,',',a,a,') <c>=',3f8.3,
     -  ' RMSD(c)=',f6.2,/,' SD(Rij)',a,':',/,(10f8.4))
1002  format('( Range of the SD scale:',f8.5,'  -',f8.5,') show')
1003  format('( Atom-atom distance SD',a,'of trajectory file ',a,
     -  ') show')
1005  format(2i5,' m')
1006  format(' NOTE: trajectory scan did not use PBC - trajectory ',
     -  'has to be centered first')
1011  format(' Range of ',a,a,' values: [',f12.5,',',f12.5,'] A')
1017  format(' SD(rij) ',a,' is zero - ',a,' items are identical')
      end
