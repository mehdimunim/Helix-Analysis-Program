      subroutine multihelix(iw0,nhx,nhxres,radtodeg,c,icaahx,
     -  icbahx,icbreshx,maxat,maxhx,mxnhx)
      dimension c(3,maxat),icaahx(maxhx,mxnhx),
     -  icbahx(mxnhx),icbreshx(mxnhx)
      parameter (MAXFRAMES=50000,MAXCOPY=600)
      common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
     -  increment,increment2,res(2,MAXFRAMES,MAXCOPY),
     -  x0(MAXCOPY),y0(MAXCOPY),nxselres,ixselres(MAXCOPY)
      parameter (MAXNHX=12,MAXNHX2=(MAXNHX*(MAXNHX-1)))
      character*2 ap_pa,in_ex
      common /signs/ itmem,normhx,isg2(MAXNHX2),memdir(MAXNHX),ap_pa(3),
     -  in_ex(2)
c     common /hxtrack/ ihxt1,ihxt2,torsav1(3,4),torsav2(3,4),x1(3),x2(3)
      dimension ax1(3),ax2(3),cent1(3),cent2(3),c12(3),ct1(3),ct2(3),
     -  ct12(3),torsats(3,4),s1(3),s2(3),e1(3),e2(3),exs1(3),
     -  exs2(3),exe1(3),exe2(3)
      real*8 dargin,dax1(3),dax2(3),dc12(3),dax1_dax2,dc12_1,dc12_2,
     -  dscprod
      character*80 line
      ixres=nhx*nhxres+1
c     write (iw0,*) 'MULTIHELIX nhxres=',nhxres
      nframes=max0(1,nframe)
      nhx2tot=(nhx*(nhx-1))/2
      i0line=0
      call blankout(line,1,80)
      do ihx=1,nhx
        incr1=(ihx-1)*nhxres
        ax1(1)=res(1,nframes,incr1+13)
        ax1(2)=res(2,nframes,incr1+13)
        ax1(3)=res(1,nframes,incr1+14)
        cent1(1)=res(2,nframes,incr1+14)
        cent1(2)=res(1,nframes,incr1+15)
        cent1(3)=res(2,nframes,incr1+15)
        halflen1=res(1,nframes,incr1+7)/2.0
        do jhx=ihx+1,nhx
          incr2=(jhx-1)*nhxres
          ax2(1)=res(1,nframes,incr2+13)
          ax2(2)=res(2,nframes,incr2+13)
          ax2(3)=res(1,nframes,incr2+14)
          cent2(1)=res(2,nframes,incr2+14)
          cent2(2)=res(1,nframes,incr2+15)
          cent2(3)=res(2,nframes,incr2+15)
          halflen2=res(1,nframes,incr2+7)/2.0
          ax1_ax2=scprod(ax1,ax2)
          ang=dacoscheck(dargin,ax1_ax2,0,6,'MULTIHELIX')*radtodeg
          dist_cc=sqrt(dist2(cent1,cent2))
          call arrdiff(cent1,cent2,c12,3)
          c12_1=scprod(c12,ax1)
          c12_2=scprod(c12,ax2)
          parfact=1.0
          if (abs(ax1_ax2) .lt. 0.01) then
c           Axes perpendicular
            a1=-c12_1
            a2=c12_2
            write (iw0,1001) ihx,jhx,'perpendicular'
          else if (abs(ax1_ax2) .gt. 0.99) then
            if (abs(ax1_ax2) .gt. 0.999) then
c             Axes parallel
              a2=0.0
              call project(cent2,cent1,ax1,ct1,a1)
              write (iw0,1001) ihx,jhx,'parallel'
              parfact=0.1
            else
              do k=1,3
                dax1(k)=ax1(k)
                dax2(k)=ax2(k)
                dc12(k)=c12(k)
              end do
              dax1_dax2=dscprod(dax1,dax2)
              dc12_1=dscprod(dc12,dax1)
              dc12_2=dscprod(dc12,dax2)
              a1=(dc12_2-dc12_1/dax1_dax2)/(1.d0/dax1_dax2-dax1_dax2)
              a2=(dc12_1-dc12_2/dax1_dax2)/(dax1_dax2-1.d0/dax1_dax2)
            end if
          else
            a1=(c12_2-c12_1/ax1_ax2)/(1.0/ax1_ax2-ax1_ax2)
            a2=(c12_1-c12_2/ax1_ax2)/(ax1_ax2-1.0/ax1_ax2)
          end if
          call paramx(cent1,ax1,a1,ct1)
          call paramx(cent2,ax2,a2,ct2)
          call arrdiff(ct1,ct2,ct12,3)
          c12norm=sqrt(scprod(ct12,ct12))
          tst1=scprod(ct12,ax1)/c12norm
          tst2=scprod(ct12,ax2)/c12norm
          if ((tst1+tst2)*parfact .gt. 0.01) then
            if (tst1+tst2 .lt. 0.1) then
              write (6,1002) 'Probable',ihx,jhx,tst1,tst2
              write (iw0,1002) 'Probable',ihx,jhx,tst1,tst2
            else
              write (6,1002) '',ihx,jhx,tst1,tst2
              write (iw0,1002) '',ihx,jhx,tst1,tst2
            end if
          end if
c         write (iw0,8766) ihx,jhx,halflen1,halflen2,a1,a2
c8766     format(' I,JHX=',2i2,' Hlen1,2=',2f6.2,' a1,2=',2f8.2)
          if (a1 .lt. -halflen1 .or. a1 .gt. halflen1 .or.
     -        a2 .lt. -halflen2 .or. a2 .gt. halflen2) then
c           Make sure the closest approach points are inside the helix
            noproj=0
            d2min=99999.0
            call paramx(cent2,ax2,-halflen2,s2)
            call paramx(cent2,ax2,halflen2,e2)
            call paramx(cent1,ax1,-halflen1,s1)
            call paramx(cent1,ax1,halflen1,e1)
            call project(s1,cent2,ax2,exs1,as1)
            call project(e1,cent2,ax2,exe1,ae1)
            call project(s2,cent1,ax1,exs2,as2)
            call project(e2,cent1,ax1,exe2,ae2)
            if (a1 .lt. -halflen1) then
c             call compared2(s1,exs1,d2min,halflen2,-halflen1,a1,a2,
c    -          as1)
              if (amax1(as2,ae2) .lt. -halflen1) noproj=1
            else if (a1 .gt. halflen1) then
c             call compared2(e1,exe1,d2min,halflen2,halflen1,a1,a2,ae1)
c             if (amin1(as1,ae1) .gt. halflen1) noproj=1
            end if
            if (a2 .lt. -halflen2) then
c             call compared2(s2,exs2,d2min,halflen1,-halflen2,a2,a1,
c    -          as2)
              if (amax1(as1,ae1) .lt. -halflen2) noproj=1
            else if (a2 .gt. halflen2) then
c             call compared2(e2,exe2,d2min,halflen1,halflen2,a2,a1,ae2)
c             if (amin1(as1,ae1) .gt. halflen2) noproj=1
            end if
            if (d2min .lt. 99999.0) then
              write (iw0,1003) ihx,jhx
            else
              if (noproj .eq. 0)
     -          write (iw0,1004) 'mismatch between',ihx,jhx
              if (noproj .eq. 1)
     -          write (iw0,1004) 'no overlap between',ihx,jhx
            end if
          end if
cn        call arrdiff(ct1,ct2,ct12,3)
cn        ct12norm=sqrt(scprod(ct12,ct12))
cn        do k=1,3
cn          ct12(k)=ct12(k)/ct12norm
cn        end do
          angrot1=0.0
          if (icbahx(ihx) .gt. 0) then
            call trnsfr(torsats(1,4),c(1,icbahx(ihx)),3)
            call project(c(1,icbahx(ihx)),cent1,ax1,torsats(1,3),a)
            call arrsum(torsats(1,3),ax1,torsats(1,2),3)
cn          call arrdiff(torsats(1,2),ct12,torsats(1,1),3)
            call project(cent2,cent1,ax1,s1,ac1)
            call arrdiff(cent2,s1,ct12,3)
            ct12norm=sqrt(scprod(ct12,ct12))
            do k=1,3
              ct12(k)=ct12(k)/ct12norm
            end do
            call arrsum(torsats(1,2),ct12,torsats(1,1),3)
            angrot1=dihangl(torsats,1,2,3,4,1,4)*radtodeg
c           write (iw0,3000)ihx,jhx,1,angrot1,((torsats(k,ii),k=1,3),
c    -        ii=1,4)
c           call trnsfr(torsav1,torsats,12)
          end if
          angrot2=0.0
          if (icbahx(jhx) .gt. 0) then
            call trnsfr(torsats(1,4),c(1,icbahx(jhx)),3)
            call project(c(1,icbahx(jhx)),cent2,ax2,torsats(1,3),a)
            call arrsum(torsats(1,3),ax2,torsats(1,2),3)
cn          call arrsum(torsats(1,2),ct12,torsats(1,1),3)
            call project(cent1,cent2,ax2,s2,ac1)
            call arrdiff(cent1,s2,ct12,3)
            ct12norm=sqrt(scprod(ct12,ct12))
            do k=1,3
              ct12(k)=ct12(k)/ct12norm
            end do
            call arrsum(torsats(1,2),ct12,torsats(1,1),3)
            angrot2=dihangl(torsats,1,2,3,4,1,4)*radtodeg
c           write (iw0,3000)ihx,jhx,2,angrot2,((torsats(k,ii),k=1,3),
c    -        ii=1,4)
c3000       format(' ihx,jhx=',2i3,' angrot',i1,'=',f8.1,' torsats=',
c    -        4(2x,3f8.2))
c           call trnsfr(torsav2,torsats,12)
          end if
c         Helix distances
          dist=sqrt(dist2(ct1,ct2))
          idist=int(50.0*dist)
          idist_cc=int(50.0*dist_cc)
          res(1,nframes,ixres)=10000*idist+idist_cc
c         End-end distances
          call dee(nframes,ixres-nhx*nhxres,ax1,ax2,cent1,cent2,
     -      halflen1,halflen2,distss,disthh,isg2ij)
          idist_ss=int(50.0*distss)
          idist_hh=int(50.0*disthh)
          res(2,nframes,ixres)=10000*idist_ss+idist_hh
c         Torsion angle over the closest aproach 'bond'
          call trnsfr(torsats(1,2),ct1,3)
          call trnsfr(torsats(1,3),ct2,3)
          call paramx(ct1,ax1,halflen1,torsats(1,1))
          call paramx(ct2,ax2,halflen2,torsats(1,4))
          dhang=dihangl(torsats,1,2,3,4,1,4)*radtodeg
c         write (77,6001) ihx,jhx
c         write (77,6002) (i,' TA ',1,(torsats(k,i),k=1,3),i=1,4)
c         write (77,6003) 'dhang',dhang
          iang=int(10.0*ang)
          idhang=int(10.0*dhang)
          res(1,nframes,nhx2tot+ixres)=10000*iang+idhang
c         Torsion angle over the center-center 'bond'
          call trnsfr(torsats(1,2),cent1,3)
          call trnsfr(torsats(1,3),cent2,3)
          call paramx(cent1,ax1,halflen1,torsats(1,1))
          call paramx(cent2,ax2,halflen2,torsats(1,4))
          dhang_cc=dihangl(torsats,1,2,3,4,1,4)*radtodeg
c         write (77,6001) ihx,jhx
c         write (77,6002) (i,' TA ',1,(torsats(k,i),k=1,3),i=1,4)
c         write (77,6003) 'dhang_cc',dhang_cc
          iamin1=icaahx(icbreshx(ihx),ihx)
          iamin2=icaahx(icbreshx(jhx),jhx)
c         if (ihx .eq. ihxt1 .and. jhx .eq. ihxt2) then
c           if (nframes .eq. 1) write (78,6001) ihx,jhx
c           write (78,6004) nframes
c           call paramx(cent1,ax1,-halflen1,x1)
c           call paramx(cent1,ax1,+halflen1,x2)
c           write (78,6002) 1,' OB1',1,x1
c           write (78,6002) 2,' OE1',1,x2
c           write (78,6002) 3,' HC1',1,cent1
c           write (78,6002) 4,' CA1',1,(c(k,iamin1),k=1,3)
c           write (78,6002) 5,' CB1',1,(c(k,icbahx(ihx)),k=1,3)
c           write (78,6002) 6,' N31',1,(torsav1(k,3),k=1,3)
c           write (78,6002) 7,' N21',1,(torsav1(k,2),k=1,3)
c           write (78,6002) 8,' N11',1,(torsav1(k,1),k=1,3)
c           write (78,6002) 9,'CL1 ',1,ct1
c           call paramx(cent2,ax2,-halflen2,x1)
c           call paramx(cent2,ax2,+halflen2,x2)
c           write (78,6002) 10,' OB2',2,x1
c           write (78,6002) 11,' OE2',2,x2
c           write (78,6002) 12,' HC2',2,cent2
c           write (78,6002) 13,' CA2',2,(c(k,iamin2),k=1,3)
c           write (78,6002) 14,' CB2',2,(c(k,icbahx(jhx)),k=1,3)
c           write (78,6002) 15,' N32',2,(torsav2(k,3),k=1,3)
c           write (78,6002) 16,' N22',2,(torsav2(k,2),k=1,3)
c           write (78,6002) 17,' N12',2,(torsav2(k,1),k=1,3)
c           write (78,6002) 18,'CL2 ',2,ct2
c           write (78,6005) 1,2
c           write (78,6005) 4,5
c           write (78,6005) 5,6
c           write (78,6005) 6,7
c           write (78,6005) 7,8
c           write (78,6005) 10,11
c           write (78,6005) 13,14
c           write (78,6005) 14,15
c           write (78,6005) 15,16
c           write (78,6005) 16,17
c           write (78,6005) 9,18
c           write (78,1000) 'ENDMDL'
c         end if
          res(2,nframes,nhx2tot+ixres)=dhang_cc
          idist_ar1=int(10.0*angrot1)
          idist_ar2=int(10.0*angrot2)
          res(1,nframes,2*nhx2tot+ixres)=10000*idist_ar1+idist_ar2
          write (iw0,1005) ihx,jhx,dist,dist_cc,distss,disthh,
     -      ap_pa(isg2ij+2)
          write (iw0,1006) ihx,jhx,ang,dhang,dhang_cc,iamin1,iamin2
          write (iw0,1007) ihx,jhx,angrot1,jhx,ihx,angrot2
          ixres=ixres+1
        end do
      end do
      return
c1000 format(a)
1001  format(' Helices',i3,' and',i3,' are ',a)
1002  format(1x,a,' PROGRAM ERROR: ihx,jhx=',2i3,' c12.m1=',f8.3,
     -  ' c12.m2=',f8.5)
1003  format(' NOTE: closest HX',i2,' - HX',i2,' contact is outside ',
     -  'the helices')
1004  format(' NOTE: ',a,' HX',i2,' - HX',i2,
     -  ' projections')
1005  format(' HXs#',i2,' - ',i3,' dist=',f7.2,' cc_dist=',f7.2,
     -  ' dee1=',f7.2,' dee2=',f7.2,1x,a)
1006  format(' HXs#',i2,' - ',i3,' ang=',f7.2,' dhang=',f7.2,
     -   ' dhang_cc=',f7.2,' CAs:',i5,i6)
1007  format(' HX',i2,' rotation wrt HX',i2,'=',f7.2,
     -  ' HX',i2,' rotation wrt HX',i2,'=',f7.2)
c6001 format('REMARK ihx,jhx=',2i3)
c6002 format('ATOM  ',i5,1x,a4,1x,'TST',1x,'D',i4,1x,3x,3f8.3,
c    -  '  1.0   0.0')
c6003 format('REMARK ',a,'=',f10.4)
c6004 format('MODEL',i6)
c6005 format('CONECT',2i5)
      end
