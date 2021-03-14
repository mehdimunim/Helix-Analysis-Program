      subroutine findbestrep(iw0,icl,iclsub,nframe,mask,irdavmn,irdmxmn,
     -  ieng,emin,emax,eav,diam,diammin,diammax,irepeng,label,llabel,
     -  iout_term,mx2d)
      dimension mask(mx2d)
      character*(*) label
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (4*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
     -  ihbpair(2,MAXBONDS),a1(MAXBONDS),a2(MAXBONDS),fill(IFILL2)
c     print *,'FINDBESTREP iw0,icl,nframe=',iw0,icl,nframe
      rdsmn=1.0e+20
      rdmn=1.0e+20
      nmem=0
      rmsdmax=0.0
      nframeused=0
      irdmxmn=0
      iravxmn=0
      do ix=1,nframe
        if (mask(ix) .gt. 0) nframeused=nframeused+1
      end do
      nsing=0
      if (icl .gt. 0) then
        if (iclsub .eq. 0) then
          if (iout_term .gt. 0)
     -      write (iout_term,1005) icl,nframeused
          if (iw0 .gt. 0) write (iw0,1005) icl,nframeused
        else
          if (iout_term .gt. 0) write (iout_term,1002)
     -      icl,iclsub,nframeused
          if (iw0 .gt. 0) write (iw0,1002) icl,iclsub,nframeused
        end if
      else
        write (6,*)
      end if
      if (ieng .gt. 0) then
        write (6,1006) icl,eav,emin,emax
        write (iw0,1006) icl,eav,emin,emax
      end if
      do ix=1,nframe
        nmem0=nmem
        if (mask(ix) .gt. 0) then
          nmem=nmem+1
          a1(ix)=0.0
          a2(ix)=0.0
          do iy=1,ix-1
            if (mask(iy) .gt. 0) then
              a1(ix)=a1(ix)+rmsd2d(iy,ix)**2
              if (rmsd2d(iy,ix) .gt. a2(ix)) a2(ix)=rmsd2d(iy,ix)
              if (rmsd2d(iy,ix) .ge. rmsdmax) then
                rmsdmax=rmsd2d(iy,ix)
                ix1=ix
                ix2=iy
              end if
            end if
          end do
          do iy=ix+1,nframe
            if (mask(iy) .gt. 0) then
              a1(ix)=a1(ix)+rmsd2d(ix,iy)**2
              if (rmsd2d(ix,iy) .gt. a2(ix)) a2(ix)=rmsd2d(ix,iy)
            end if
          end do
          if (a1(ix) .lt. rdsmn) then
            rdsmn=a1(ix)
            irdavmn=ix
          end if
          if (a2(ix) .lt. rdmn) then
            rdmn=a2(ix)
            irdmxmn=ix
          end if
        end  if
        if (nmem .gt. nmem0) a1(ix)=a1(ix)/float(nmem-nmem0)
      end do
      if (nmem .gt. 1) then
        if (iout_term .gt. 0)
     -    write (iout_term,1010) irdmxmn,label(1:llabel),rdmn,
     -      label(1:llabel),rmsdmax,ix1,ix2,label(1:llabel),
     -      label(1:llabel),ix1,irdmxmn,rmsd2d(ix1,irdmxmn),
     -      label(1:llabel),ix2,irdmxmn,rmsd2d(ix2,irdmxmn)
        if (iw0 .gt. 0)
     -    write (iw0,1010) irdmxmn,label(1:llabel),rdmn,
     -    label(1:llabel),rmsdmax,ix1,ix2,label(1:llabel),
     -    label(1:llabel),ix1,irdmxmn,rmsd2d(ix1,irdmxmn),
     -    label(1:llabel),ix2,irdmxmn,rmsd2d(ix2,irdmxmn)
        if (label(1:4) .eq. 'RMSD') then
          if (iout_term .gt. 0)
     -      write (iout_term,1011) irdavmn,rdsmn/nframeused,ix1,irdavmn,
     -      rmsd2d(ix1,irdavmn),ix2,irdavmn,rmsd2d(ix2,irdavmn)
          if (iw0 .gt. 0)
     -      write (iw0,1011) irdavmn,rdsmn/nframeused,ix1,irdavmn,
     -      rmsd2d(ix1,irdavmn),ix2,irdavmn,rmsd2d(ix2,irdavmn)
          if (diammin .gt. rmsdmax) diammin=rmsdmax
          if (diammax .lt. rmsdmax) diammax=rmsdmax
          diam=rmsdmax
        end if
        if (irdavmn .ne. irdmxmn) then
          if (iout_term .gt. 0)
     -      write (iout_term,1001) label(1:llabel),irdavmn,irdmxmn,
     -        rmsd2d(irdavmn,irdmxmn)
          write (iw0,1001) label(1:llabel),irdavmn,irdmxmn,
     -      rmsd2d(irdavmn,irdmxmn)
        end if
      end if
      if (label(1:4) .eq. 'RMSD' .and. iw0 .gt. 0) then
        write (iw0,1003)
        do ix=1,nframe
          if (mask(ix) .gt. 0) then
            write (iw0,1004) ix,sqrt(a1(ix)),a2(ix)
          end if
        end do
      end if
      return
1001  format(1x,a,' between the two center estimates',
     -  ' (',i5,',',i5,')=',f8.2)
1002  format(/,' Cluster #',i4,' Subcluster #',i3,
     -  ' Number of members=',i5)
1003  format(/,' RMSD average and maximum (over other cluster members)',
     -  ' for each frame analyzed')
1004  format(' Frame #',i4,' sqrt(<MSD>)=',f12.3,' Maximum RMSD=',f8.2)
1005  format(/,' Cluster #',i4,' Number of members=',i5)
1006  format(' Cluster #',i4,' Average energy: ',e12.5,
     -  ' Range: [',e12.5,',',e12.5,']')
1011  format(' Cluster center based on the lowest mean MSD:    #',i5,
     -  ' (<MSD>=',f8.2,')',/
     -  ' Cluster radius based on the lowest mean MSD:   ',/,
     -  16x,'RMSD(',i5,',',i5,')=',f8.2,'  RMSD(',i5,',',i5,')=',f8.2)
1010  format(' Cluster center: #',i5,' - based on the lowest maximum ',
     -  a,' (',f8.2,')',/
     -  ' Largest ',a,' (cluster diameter) is ',f8.2,
     -  ' between #',i5,' & #',i5,/,
     -  ' Cluster radius based on the lowest maximum ',a,':',/,
     -  1x,a,'(',i5,',',i5,')=',f8.2,2x,a,'(',i5,',',i5,')=',f8.2)
      end
