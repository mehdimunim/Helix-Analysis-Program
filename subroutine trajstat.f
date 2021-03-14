      subroutine trajstat(iout,ndials,maxdials,angname,langname,nrdat,
     -  maxrdat,irescolf,rname,lrname,corr12,incrres,icorr,now6,
     -  title,ltitle,radtodeg)
      character*(*) angname(maxdials),rname(maxrdat),title
      dimension langname(maxdials),lrname(maxrdat)
      parameter (MAXFRAMES=50000,MAXCOPY=600)
      common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
     -  increment,increment2,res(2,MAXFRAMES,MAXCOPY),
     -  x0(MAXCOPY),y0(MAXCOPY),nxselres,ixselres(MAXCOPY)
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-MAX2D*MAX2D)
      common /nnwork/ ccc(MAX2D,MAX2D),fill(IFILL2)
      dimension cv(MAXCOPY),av_r(MAXCOPY),sd_r(MAXCOPY),lcorrtyp(2)
      real*8 ddd,sum,sum2,sinsum,cossum,sinsum_ab,sinsum2,cossum2,
     -  sinsum2_a,sinsum2_b,sin_ai,sin_bi,av_a(MAXCOPY),ang_ai,ang_aj,
     -  ang_bi,ang_bj,sin_aij,r12,r13,r23,csd,ssd,csc,crs,crc
      character*26 corrtyp(2)
      data corrtyp /'Fisher & Lee','Jammalamadaka and SenGupta'/,
     -  lcorrtyp /12,26/
      write (iout,*)
      if (ltitle .gt. 1) write (iout,*) title(1:ltitle)
      if (maxdials .gt. MAXCOPY-1) then
        write (6,2004) maxdials,MAXCOPY
        stop
      end if
      if (nframe .lt. 1) return
      do id=1,ndials
        sinsum=0.d0
        cossum=0.d0
        do i=1,nframe
          cossum=cossum+res(1,i,incrres+id)
          sinsum=sinsum+res(2,i,incrres+id)
        end do
        sum=dsqrt(sinsum**2+cossum**2)
        cv(id)=1.0-sum/nframe
        cossum=cossum/sum
        av_a(id)=dacoscheck(cossum,csum,1,6,'TRAJSTAT')
        if (sinsum .lt. 0.d0) av_a(id)=-av_a(id)
        if (now6 .eq. 0) write (6,2078) angname(id)(1:langname(id)),
     -    av_a(id)*radtodeg,cv(id)
        write (iout,2078) angname(id)(1:langname(id)),
     -    av_a(id)*radtodeg,cv(id)
      end do
      do irdat=1,nrdat
        id=irescolf+(irdat-1)/2
        ic=mod(irdat-1,2)+1
        sum=0.d0
        sum2=0.d0
        do i=1,nframe
          sum=sum+res(ic,i,incrres+id)
          sum2=sum2+res(ic,i,incrres+id)**2
        end do
        av_r(irdat)=sum/nframe
        sd_r(irdat)=dsqrt(dabs(sum2/nframe-av_r(irdat)**2))
        if (now6 .eq. 0) write (6,2077) rname(irdat)(1:lrname(irdat)),
     -    av_r(irdat),sd_r(irdat)
        write (iout,2077)
     -    rname(irdat)(1:lrname(irdat)),av_r(irdat),sd_r(irdat)
      end do
      if (icorr .gt. 0) then
        if (ndials .gt. 0) then
c         Calculate circular correlation coefficient between angles
          if (ndials .le. 10 .and. now6 .eq. 0)
     -      write (6,2001) corrtyp(icorr)(1:lcorrtyp(icorr))
          write (iout,2001) corrtyp(icorr)(1:lcorrtyp(icorr))
          do id=1,ndials
            ccc(id,id)=1.0
            do jd=id+1,ndials
              sinsum_ab=0.d0
              sinsum2_a=0.d0
              sinsum2_b=0.d0
              if (icorr .eq. 1) then
c               Fisher & Lee (slow O(nframe^2))
                do i=1,nframe
                  ang_ai=dacoscheck(ddd,res(1,i,incrres+id),0,6,
     -              'TRAJSTAT')
                  if (res(2,i,incrres+id) .lt. 0.0) ang_ai=-ang_ai
                  ang_bi=dacoscheck(ddd,res(1,i,incrres+jd),0,6,
     -              'TRAJSTAT')
                  if (res(2,i,incrres+jd) .lt. 0.0) ang_bi=-ang_bi
                  do j=i+1,nframe
                    ang_aj=dacoscheck(ddd,res(1,j,incrres+id),0,6,
     -                'TRAJSTAT')
                    if (res(2,j,incrres+id) .lt. 0.0) ang_aj=-ang_aj
                    ang_bj=dacoscheck(ddd,res(1,j,incrres+jd),0,6,
     -                'TRAJSTAT')
                    if (res(2,j,incrres+jd) .lt. 0.0) ang_bj=-ang_bj
                    sin_aij=dsin(ang_ai-ang_aj)
                    sin_bij=dsin(ang_bi-ang_bj)
                    sinsum_ab=sinsum_ab+sin_aij*sin_bij
                    sinsum2_a=sinsum2_a+sin_aij**2
                    sinsum2_b=sinsum2_b+sin_bij**2
                  end do
                end do
              else
c               Jammalamadaka and SenGupta (2001) (fast, O(nframe))
                do i=1,nframe
                  ang_ai=dacoscheck(ddd,res(1,i,incrres+id),0,6,
     -              'TRAJSTAT')
                  if (res(2,i,incrres+id) .lt. 0.0) ang_ai=-ang_ai
                  ang_bi=dacoscheck(ddd,res(1,i,incrres+jd),0,6,
     -              'TRAJSTAT')
                  if (res(2,i,incrres+jd) .lt. 0.0) ang_bi=-ang_bi
                  sin_ai=dsin(ang_ai-av_a(id))
                  sin_bi=dsin(ang_bi-av_a(jd))
                  sin_aij=sin_ai*sin_bi
                  sinsum_ab=sinsum_ab+sin_aij
                  sinsum2_a=sinsum2_a+sin_ai**2
                  sinsum2_b=sinsum2_b+sin_bi**2
                end do
              end if
              ccorr=(sinsum_ab)/dsqrt(sinsum2_a*sinsum2_b)
              if (ndials .le. 10 .and. now6 .eq. 0)
     -          write (6,2000) angname(id)(1:langname(id)),
     -            angname(jd)(1:langname(jd)),ccorr
              write (iout,2000) angname(id)(1:langname(id)),
     -          angname(jd)(1:langname(jd)),ccorr
              ccc(id,jd)=ccorr
              ccc(jd,id)=ccorr
            end do
          end do
        end if
        if (nrdat .gt. 0) then
c         Calculate standard correlation coefficient between regular data
          if (now6 .eq. 0) write (6,2002) ' '
          write (iout,2002) ' '
          do irdat=1,nrdat
            id=irescolf+(irdat-1)/2
            ic=mod(irdat-1,2)+1
            do jrdat=irdat+1,nrdat
              jd=irescolf+(jrdat-1)/2
              jc=mod(jrdat-1,2)+1
              call correl(iout,res(1,1,incrres+id),ic,
     -          rname(irdat),lrname(irdat),av_r(irdat),sd_r(irdat),
     -          res(1,1,incrres+jd),jc,rname(jrdat),lrname(jrdat),
     -          av_r(jrdat),sd_r(jrdat),corr,nframe,0,now6)
              if (irdat .eq. 1 .and. jrdat .eq. 2) corr12=corr
            end do
          end do
        end if
        if (ndials*nrdat .gt. 0) then
c         Calculate linear-circular correlation coefficient between
c         angles and regular data (Mardia)
          if (ndials .le. 10 .and. now6 .eq. 0) write (6,2003)
          write (iout,2003)
          do id=1,ndials
            do jrdat=1,nrdat
              corrsum=0.d0
              jd=irescolf+(jrdat-1)/2
              jc=mod(jrdat-1,2)+1
              sinsum=0.d0
              sinsum2=0.d0
              cossum=0.d0
              cossum2=0.d0
              crc=0.d0
              crs=0.d0
              csc=0.d0
              do i=1,nframe
                ang_ai=dacoscheck(ddd,res(1,i,incrres+id),0,6,
     -            'TRAJSTAT')
                if (res(2,i,incrres+id) .lt. 0.0) ang_ai=-ang_ai
                c=res(1,i,incrres+id)
                s=res(2,i,incrres+id)
                x=res(jc,i,incrres+jd)
                crc=crc+x*c
                crs=crs+x*s
                csc=csc+s*c
                sinsum=sinsum+s
                cossum=cossum+c
                sinsum2=sinsum2+s**2
                cossum2=cossum2+c**2
                corrsum=corrsum+ang_ai*x
              end do
              corrsum=corrsum*radtodeg
              cav=cossum/nframe
              sav=sinsum/nframe
              csd=(cossum2/nframe-cav**2)
              ssd=(sinsum2/nframe-sav**2)
              if (csd .gt. 0.d0 .and. ssd .gt. 0.d0) then
                csd=dsqrt(csd)
                ssd=dsqrt(ssd)
                r12=(crc/nframe-cav*av_r(jrdat))/(csd*sd_r(jrdat))
                r13=(crs/nframe-sav*av_r(jrdat))/(ssd*sd_r(jrdat))
                r23=(csc/nframe-sav*cav)/(csd*ssd)
                ccorr=(r12**2+r13**2-2.0*r12*r13*r23)/(1-r23**2)
              else
                ccorr=0.0
              end if
              if (ndials .le. 10 .and. now6 .eq. 0)
     -          write (6,2000) angname(id)(1:langname(id)),
     -            rname(jrdat)(1:lrname(jrdat)),ccorr
              write (iout,2000) angname(id)(1:langname(id)),
     -          rname(jrdat)(1:lrname(jrdat)),ccorr
            end do
          end do
        end if
      end if
c     write (77,*) 'ndials=',ndials
c     do i=1,ndials
c       write (77,1004) (ccc(i,j),j=1,ndials)
c     end do
c1004 format(5e13.6)
      return
2000  format(1x,a,' - ',a25,' ccc=',f8.5)
2001  format(' ccc: Circular correlation coefficient (',a,')')
2002  format(/,a,'Pearson correlation coefficient')
2003  format(' ccc: Linear-circular correlation coefficient (Mardia)')
2004  format(' PROGRAM ERROR: argument maxdials of subroutine trajstat',
     -  ' exceeds',/,' the parameter MAXCOPY (',i4,' vs ',i4,')')
2077  format(1x,a,' average=',f10.2,' S.D.=',f8.2,a,'CV=',f8.5)
2078  format(1x,a,' average=',f10.2,' CV=',f8.5)
      end
