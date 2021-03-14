      subroutine reportclust(ndim,icl0,nfclst,nlclst,ifclst,ilclst,
     -  index2d,value,it1,ifa_s,ila_s,ih,cv,indexa,irepav,irepmx,
     -  irepeng,irepkm,engcl,nhbdist,etotsaved,ietotsaved,ifindbestrep,
     -  label,llabel,isorttype,idistprint,nomemprint,iout,maxframe,mx2d)
      dimension it1(mx2d),ifclst(mx2d),ilclst(mx2d),index2d(mx2d),
     -  indexa(mx2d),value(mx2d),ifa_s(mx2d),ila_s(mx2d),ih(mx2d),
     -  irepav(mx2d),irepmx(mx2d),irepeng(mx2d),irepkm(mx2d),
     -  nhbdist(mx2d),cv(mx2d),engcl(mx2d),etotsaved(2,maxframe)
      character*(*) label
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (6*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
     -  ihbtores(MAXBONDS),nusepair(MAXBONDS),nhb_atot(MAXBONDS),
     -  nhb_rtot(MAXBONDS),ixrank(MAXBONDS),rankav(MAXBONDS),
     -  fill(IFILL2)
      dimension diamcl(MAX2D)
      data emax /0.0/
c     print *,'REPORTCLUST nfclst,nlclst,isorttype,ifindbestrep=',
c    -  nfclst,nlclst,isorttype,ifindbestrep
c      write (77,8822) (index2d(i),i=1,ilclst(nlclst))
c8822  format(' index2d:',25i5)
c     do i=1,ndim
c       write (77,9782) i,(rmsd2d(i,j),j=1,ndim)
c9782   format(i4,' ccc=',30f5.2)
c     end do
      nclust=nlclst-nfclst+1
      do ic=nfclst,nlclst
        ixsum=0
        do i=ifclst(ic),ilclst(ic)
          ixsum=ixsum+index2d(i)
        end do
        nmem=ilclst(ic)-ifclst(ic)+1
        rankav(ic)=float(ixsum)/float(nmem)
      end do
      call indexit(ixrank,1,nclust,0)
      call mrgsrt(6,ixrank,rankav(nfclst),nclust,ifa_s,ila_s,ih,cv,
     -  max2d)
      do icc=1,nclust
        ifa_s(icc)=ifclst(ixrank(icc)+nfclst-1)
        ila_s(icc)=ilclst(ixrank(icc)+nfclst-1)
      end do
      do ic=nfclst,nlclst
        ifclst(ic)=ifa_s(ic-nfclst+1)
        ilclst(ic)=ila_s(ic-nfclst+1)
      end do
      iout_term=6
      if (nlclst-nfclst .gt. 9) then
        iout_term=0
        print *,'Cluster characteristics are printed on the .rd2 file ',
     -    'only'
      end if
      cdiammin=100000.0
      cdiammax=0.0
      if (ifindbestrep .gt. 0) then
        call indexit(indexa,1,ndim,0)
        write (6,1015)
        call findbestrep(0,0,0,ndim,indexa,irdavmn,irdmxmn,
     -    0,0.0,0.0,0.0,cdiam,cdiammin,cdiammax,0,'RMSD',4,0,MAX2D)
      end if
      maxmem=0
      diammin=100000.0
      diammax=0.0
      do ic=nfclst,nlclst
        nmem=ilclst(ic)-ifclst(ic)+1
        if (nmem .gt. maxmem) maxmem=nmem
        if (nomemprint .eq. 0)
     -    write (iout,2090) ic,nmem,rankav(ic)
        if (isorttype .eq. 1) then
c         Sort by increasing member index
          do i=ifclst(ic),ilclst(ic)
            value(i-ifclst(ic)+1)=index2d(i)
          end do
          call indexit(it1,1,nmem,0)
          call mrgsrt(6,it1,value,nmem,ifa_s,ila_s,ih,cv,max2d)
          do i=ifclst(ic),ilclst(ic)
            index2d(i)=value(i-ifclst(ic)+1)
          end do
          i00=ifclst(ic)-1
          i0=0
          do while (i0 .lt. nmem)
            if (nomemprint .eq. 0) write (iout,2091)ic,(index2d(i00+i),
     -        i=i0+1,min0(i0+10,nmem))
            i0=i0+10
          end do
        else if (isorttype .gt. 1) then
c         Sort by decreasing occurrence
          do i=ifclst(ic),ilclst(ic)
            it1(i-ifclst(ic)+1)=index2d(i)
            value(i-ifclst(ic)+1)=-nhbdist(index2d(i))
          end do
          call mrgsrt(6,it1,value,nmem,ifa_s,ila_s,ih,cv,max2d)
          call trnsfi(index2d(ifclst(ic)),it1,nmem)
          i0=0
          do while (i0 .lt. nmem)
            if (nomemprint .eq. 0) write (iout,2091) ic,(it1(i),
     -        i=i0+1,min0(i0+10,nmem))
            i0=i0+10
          end do
        end if
        emax=-1.e+32
        if (ietotsaved .gt. 0) then
          emin=1.e+32
          esum=0.0
          do i=ifclst(ic),ilclst(ic)
            etot=etotsaved(1,index2d(i))
            esum=esum+etot
            if (emin .gt. etot) then
              emin=etot
              irepeng(ic)=index2d(i)
            end if
            if (emax .lt. etot) emax=etot
          end do
          engcl(ic)=esum/(ilclst(ic)-ifclst(ic)+1)
        end if
        if (ifindbestrep .eq. 1) then
          call zeroiti(indexa,0,ndim)
          do i=ifclst(ic),ilclst(ic)
            indexa(index2d(i))=1
          end do
          if (icl0 .eq. 0) then
            call findbestrep(iout,ic,0,ndim,indexa,irepav(ic),
     -        irepmx(ic),ietotsaved,emin,emax,engcl(ic),diamcl(ic),
     -        diammin,diammax,irepeng(ic),label,llabel,iout_term,max2d)
          else
            call findbestrep(iout,icl0,ic,ndim,indexa,irepav(ic),
     -        irepmx(ic),ietotsaved,emin,emax,engcl(ic),diamcl(ic),
     -        diammin,diammax,irepeng(ic),label,llabel,iout_term,max2d)
          endif
        end if
      end do
      write (6,2095) (ilclst(ic)-ifclst(ic)+1,ic=1,nclust)
      write (iout,2095) (ilclst(ic)-ifclst(ic)+1,ic=1,nclust)
      if (ifindbestrep .eq. 1) then
        write (6,2101) (diamcl(ic),ic=1,nclust)
        write (iout,2101) (diamcl(ic),ic=1,nclust)
      end if
      if (maxmem .gt. 1 .and. ifindbestrep .gt. 0) then
        write (6,2099)  cdiammin,diammin,diammax
        write (iout,2099)  cdiammin,diammin,diammax
      end if
      if (ifindbestrep .eq. 1) then
        if (label(1:4) .eq. 'RMSD') then
          write (iout,2096) 'average RMSD'
          i0=0
          do while (i0 .lt. nclust)
            write (iout,2098) ' clrepa ',
     -        (irepav(i),i=i0+1,min0(i0+10,nclust))
            i0=i0+10
          end do
        end if
        write (iout,2097) 'maximum ',label(1:llabel)
        i0=0
        do while (i0 .lt. nclust)
          write (iout,2098) ' clrepm ',
     -      (irepmx(i),i=i0+1,min0(i0+10,nclust))
          i0=i0+10
        end do
        write (iout,2089)
        i0=0
        do while (i0 .lt. nclust)
          write (iout,2098) ' clfrst ',
     -      (index2d(ifclst(i)),i=i0+1,min0(i0+10,nclust))
          i0=i0+10
        end do
        if (ietotsaved .gt. 0 .and. ifindbestrep .eq. 1) then
          write (iout,2096) 'lowest energy'
          i0=0
          do while (i0 .lt. nclust)
            write (iout,2094) ' clrepe ',
     -        (irepeng(ic),etotsaved(1,irepeng(ic)),
     -        ic=i0+1,min0(i0+3,nclust))
            i0=i0+3
          end do
          write (iout,2100) (engcl(ic),ic=1,nclust)
        end if
        nz=0
        do i=1,nclust
          if (irepkm(i) .eq. 0) nz=nz+1
        end do
        if (nz .eq. 0) then
          write (iout,2096) 'k-medoids clustering centers'
          i0=0
          do while (i0 .lt. nclust)
            write (iout,2098) ' clrepkm',
     -        (index2d(irepkm(i)),i=i0+1,min0(i0+10,nclust))
            i0=i0+10
          end do
        end if
        if (idistprint .eq. 1 .and. nlclst-nfclst .gt. 0) then
          do ic=nfclst,nlclst
            if (label(1:4) .eq. 'RMSD') then
              write (iout,2080) ic,
     -          (rmsd2d(irepav(ic),irepav(jc)),jc=nfclst,nlclst)
            else
              write (iout,2092) ic,label(1:llabel),
     -          (rmsd2d(irepmx(ic),irepmx(jc)),jc=nfclst,nlclst)
            end if
          end do
        end if
      end if
      return
1015  format(1x,a)
2080  format(' RMSDs between <MSD>-based center of cluster',i4,' and ',
     -  /,' the other cluster centers:',/,(10f8.2))
2089  format(/,' First member of each cluster:')
2090  format(/,' Cluster #',i4,' contains ',i4,' members (<rank>=',
     -  f8.1,'):')
2091  format(' clmem ',i4,':',10i6)
2092  format(' Cluster',i4,1x,a,'-based center-center distances:',/,
     -  (10f8.2))
2094  format(a,3(i8,' (',e12.5,')'))
2095  format(/,' Number of members in each cluster:',/,(10i5))
2096  format(/,' Cluster representatives (centers) based on ',a,':')
2097  format(/,' Cluster representatives (centers) based on ',a,a,':')
2098  format(1x,a,':',10i6)
2099  format(' Diameter of the whole set:',f8.2,
     -  ' Cluster diameter range: [',f8.2,',',f8.2,']')
2100  format(/,' Average energy of the cluster members:',/,(5e13.5))
2101  format(' Cluster diameters:',/,(10f8.2))
      end
