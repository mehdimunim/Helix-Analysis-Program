      subroutine clusterdistr(ndim,iout,rdlim,rdmin,rdmax,nhbdist,
     -  it1,it2,it3,it4,ifclst,ilclst,nrdclust,indexa,index2d,ixclst,
     -  it2d,value,ifa_s,ila_s,ih,cv,rdclustdef,rdclust,etotsaved,
     -  ietotsaved,label,llabel,isorttype,ifindbestrep,irepav,irepmx,
     -  irepeng,irepkm,engcl,c,cent,cent_prev,maxct,ihelp,iclstyp,
     -  iwt,nomemprint,label2d,llabel2d,ilabel2d,idistp,nosortindex2d,
     -  iuseindex2d,mx2d,maxframe)
      dimension rdlim(mx2d),nhbdist(mx2d),it1(mx2d),it2(mx2d),
     -  it3(mx2d),it4(mx2d),ifclst(mx2d),ilclst(mx2d),ixclst(mx2d),
     -  it2d(mx2d),index2d(mx2d),indexa(mx2d),value(mx2d),ifa_s(mx2d),
     -  ila_s(mx2d),ih(mx2d),irepav(mx2d),irepmx(mx2d),iwt(mx2d),
     -  cv(mx2d),etotsaved(2,maxframe),irepeng(mx2d),irepkm(mx2d),
     -  engcl(mx2d),c(3,ndim),cent(3,maxct),cent_prev(3,maxct)
      character*(*) label
      character*(*) label2d(mx2d)
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      dimension cutofflist(5),imemcut(MAX2D,5),ifcl_prev(MAX2D),
     -  ilcl_prev(MAX2D),icent_fin(MAX2D)
      character*18 memlab
      character*200 memdir,memfilename
      character*500 line
      character*1 ans
      character*41 clstyp
      common /cluster_typ/ nclstyp,inumclst(9),ireadcutoff(9),
     -  lclstyp(9),clstyp(9)
      data iu_clst /60/
c     Clustering
c     print *,'CLUSTERD ndim,rdmin,rdmax,ietotsaved=',
c    -  ndim,rdmin,rdmax,ietotsaved
c      write (iout,7934) (index2d(i),i=1,ndim)
c7934  format(' CLUSTERDISTR INDEX2D=',15i4)
      idistprint=idistp
      nrdclust=0
      irepeng(1)=0
      nnmin=0
      idenclstyp=0
      if (rdmax .le. 0.0) then
        print *,'ERROR: largest distance is nonpositive:',rdmax
        return
      end if
      iclstyp=1
8910  call quiz(ans,iclstyp,' ',' ',0,'clustering algorithm',20,
     -  0,5,6,40)
      if (ans .eq. 'q') then
        if (nrdclust .eq. 0) then
          print *,'NOTE: no clustering was done'
          call askyn('No clustering was done - still want to quit',43,
     -      1,1,iquit,130,0)
          if (iquit .eq. 0) go to 8910
        end if
        return
      end if
      if (maxct .le. 1 .and. iclstyp .eq. 5) then
        write (6,*)
     -   'COM as cluster center only applies to clustering atoms'
        go to 8910
      end if
      write (iout,2060) clstyp(iclstyp),ndim
      write (6,2060) clstyp(iclstyp),ndim
      if (iclstyp .ne. 4) call zeroiti(irepkm,0,ndim)
      if (inumclst(iclstyp) .eq. 1) then
        call getint('Number of clusters requested',28,999999,1,ndim,
     -    nrdclust,ihelp)
        rdclust=0.0
        write (iout,2061) nrdclust
      else if (ireadcutoff(iclstyp) .eq. 1) then
        write (line,1014) label(1:llabel)
        call getreal(line,22+llabel,rdclustdef,rdclust,1,ihelp)
        nhbcorrclust=1
      end if
      if (iclstyp .eq. 9) then
c       Set density type and neighbor min
        call quiz(ans,idenclstyp,' ',' ',0,'density clustering variant',
     -    26,0,5,6,141)
        if (idenclstyp .lt. 2) call getint(
     -    'Minimum number of common neighbors',34,999999,1,ndim,nnmin,
     -    000)
        if (idenclstyp .eq. 3) call getint(
     -    'Minimum number of individual neighbors',38,999999,1,ndim,
     -    nnmin,000)
      end if
      write (iout,2062)
      if (iclstyp .ne. 3 .and. idistprint .ge. 0)
     -  call askyn('Do you want to print all cluster-cluster distances',
     -    50,1, -1,idistprint,130,0)
      call zeroit(rdlim,ndim)
      rdlim(1)=rdmax
      rdlim(ndim)=rdmin
      if (rdlim(ndim) .eq. 0.0) rdlim(ndim)=0.000001
      call trnsfi(it2d,index2d,ndim)
      if (iclstyp .eq. 1 .or. iclstyp .eq. 4 .or. iclstyp .eq. 5 .or.
     -    iclstyp .eq. 6 .or. iclstyp .eq. 8 .or. iclstyp .eq. 9) then
c       Just cluster using rdclust as the threshold
        if (ans .ne. 'k')
     -    write (iout,2092) label(1:llabel),rdclust
        call rmsdcluster(rdclust,1,ndim,index2d,iwt,ixclst,ifclst,
     -    ilclst,it1,it2,it3,irepav,irepmx,it4,c,cent,cent_prev,irepkm,
     -    maxct,iclstyp,idenclstyp,nnmin,nrdclust,nosortindex2d,
     -    iuseindex2d,ifail,1,label,llabel,iout)
c      write (iout,7824) (ifclst(i),ilclst(i),i=1,nrdclust)
c7824  format(' After RMSDCLUSTER Cluster limits: ',('[',i5,',',i5,']'))
c      write (iout,6792) 'IXCLST',(ixclst(i),i=1,ndim)
c6792  format(' AFTER RMSDCLUSTER ',a,':',/,(20i5))
        if (ifail .gt. 0) go to 8910
        nclust=nrdclust
c       Members of cluster ic: (index2d(i),i=ifclst(ic),ilclst(ic))
        if (ans .ne. 'k') write (6,2090) label(1:llabel),rdclust,nclust
        call reportclust(ndim,0,1,nclust,ifclst,ilclst,index2d,value,
     -    it1,ifa_s,ila_s,ih,cv,indexa,irepav,irepmx,irepeng,irepkm,
     -    engcl,nhbdist,etotsaved,ietotsaved,ifindbestrep,label,llabel,
     -    isorttype,idistprint,nomemprint,iout,maxframe,mx2d)
      else if (iclstyp .eq. 2 .or. iclstyp .eq. 7) then
c       Vary the threshold until nrdclust clusters result
        write (6,2093) nrdclust,rdmin,rdmax
        write (iout,2093) nrdclust,rdmin,rdmax
        call zeroiti(indexa,0,ndim)
        nclust=1
        rdmn=rdmin
        rdmx=rdmax
        nclmx=nrdclust
        do while (rdlim(nclmx) .eq. 0.0)
          nclmx=nclmx-1
        end do
        rdmx=rdlim(nclmx)
        nclmn=nrdclust
        do while (rdlim(nclmn) .eq. 0.0)
          nclmn=nclmn+1
        end do
        rdmn=rdlim(nclmn)
        ntry=0
        call trnsfi(it4,index2d,ndim)
        do while (nclust .ne. nrdclust .and. ntry .le. ndim)
          call trnsfi(index2d,it4,ndim)
          rdclust=(rdmn+rdmx)/2.0
          call rmsdcluster(rdclust,1,ndim,index2d,iwt,ixclst,ifclst,
     -      ilclst,it1,it2,it3,irepav,irepmx,it4,c,cent,cent_prev,
     -      icent_fin,maxct,iclstyp,idenclstyp,nnmin,nclust,
     -      nosortindex2d,iuseindex2d,ifail,ntry,label,llabel,iout)
          if (ifail .gt. 0) go to 8910
          if (nclust .lt. nrdclust) rdmx=rdclust
          if (nclust .gt. nrdclust) rdmn=rdclust
          write (6,2090) label(1:llabel),rdclust,nclust
          write (iout,2090) label(1:llabel),rdclust,nclust
          if (indexa(nclust) .gt. 25) then
            write (iout,2094) nrdclust,nclust
            write (6,2094) nrdclust,nclust
            nrdclust=nclust
            ntry=ndim
          end if
          indexa(nclust)=indexa(nclust)+1
          rdlim(nclust)=rdclust
          ntry=ntry+1
        end do
        write (iout,2066) label(1:llabel),nclust,rdclust
        write (6,2066) label(1:llabel),nclust,rdclust
        if (nclust .ne. nrdclust) then
          write (iout,2094) nrdclust,nclust
          write (6,2094) nrdclust,nclust
          nrdclust=nclust
        end if
        call reportclust(ndim,0,1,nclust,ifclst,ilclst,index2d,value,
     -    it1,ifa_s,ila_s,ih,cv,indexa,irepav,irepmx,irepeng,irepkm,
     -    engcl,nhbdist,etotsaved,ietotsaved,ifindbestrep,label,llabel,
     -    isorttype,idistprint,nomemprint,iout,maxframe,mx2d)
      else if (iclstyp .eq. 3) then
100     call getint('Number of cutoffs',17,5,1,5,ncutoff,00)
        if (ncutoff .le. 1) go to 100
        call askyn('Do you want uniformly spaced cutoffs',36,1,1,
     -    iunif,0,0)
        if (iunif .eq. 1) then
          call transform_dist(rdmin,rdtmin)
          call transform_dist(rdmax,rdtmax)
          if (rdtmin .lt. rdtmax) then
            call getreal('Largest cutoff',14,rdtmax,rcutulim,1,00)
            call getreal('Smallest cutoff',15,rdtmin,rcutllim,1,00)
          else
            call getreal('Largest cutoff',14,rdtmin,rcutulim,1,00)
            call getreal('Smallest cutoff',15,rdtmax,rcutllim,1,00)
          end if
          do icut=1,ncutoff
            cutofflist(icut)=
     -        rcutulim-(icut-1)*(rcutulim-rcutllim)/(ncutoff-1)
          end do
        else
          print *,'Specify the cutoffs (in decreasing order)'
          rcutprev=rdmax+1.0
          do icut=1,ncutoff
110         write (line(1:11),1015) icut
            call getreal(line,11,999999.0,cutofflist(icut),1,00)
            if (cutofflist(icut) .ge. rdmax) then
              print *,'Cutoff read exceeds the largest distance'
              go to 110
            else if (cutofflist(icut) .ge. rcutprev) then
              print *,'Cutoffs need to be in decreasing order'
              go to 110
            end if
            rcutprev=cutofflist(icut)
          end do
        end if
        write (6,1016) (cutofflist(icut),icut=1,ncutoff)
        write (iout,*) 'Hierarchical clustering with multiple cutoffs'
        write (iout,1016) (cutofflist(icut),icut=1,ncutoff)
        call getint('Cluster level to calculate average',34,ncutoff,1,
     -    ncutoff,lev_avg,00)
        call askyn('Do you want to write cluster member files',41,1,-1,
     -    memfiles,129,0)
        if (memfiles .eq. 1) then
111       call getname(memdir,lmemdir,'Name of the directory to write',
     -      30,60,'',0,0,0,0)
          call checkdir(memdir,lmemdir,iu_clst,iopen)
          if (iopen .gt. 0) then
            write (6,1024) memdir(1:lmemdir)
            go to 111
          end if
          memdir(lmemdir+1:lmemdir+2)='/C'
          lmemdir=lmemdir+2
          memfilename=memdir
        end if
        write (iout,1023) lev_avg
        ifcl_prev(1)=1
        ilcl_prev(1)=ndim
        ncl_prev=1
        call indexit(index2d,1,ndim,0)
        write (iout,1019)
        do icut=1,ncutoff
          ncltot=0
          do icl=1,ncl_prev
c           write (iout,*)
c    -        'cutofflist(icut),ifcl_prev(icl),ilcl_prev(icl)=',
c    -        cutofflist(icut),ifcl_prev(icl),ilcl_prev(icl)
            call rmsdcluster(cutofflist(icut),ifcl_prev(icl),
     -        ilcl_prev(icl),index2d,iwt,ixclst,ifclst(ncltot+1),
     -        ilclst(ncltot+1),it1,it2,it3,irepav,irepmx,it4,c,cent,
     -        cent_prev,icent_fin,maxct,iclstyp,idenclstyp,nnmin,nclust,
     -        nosortindex2d,iuseindex2d,ifail,1,label,llabel,iout)
            if (ifail .gt. 0) go to 8910
            do ic=ncltot+1,ncltot+nclust
              do ia=ifclst(ic),ilclst(ic)
                imemcut(index2d(ia),icut)=ic-ncltot
              end do
            end do
            ncltot=ncltot+nclust
          end do
          call trnsfi(ifcl_prev,ifclst,ncltot)
          call trnsfi(ilcl_prev,ilclst,ncltot)
          ncl_prev=ncltot
        end do
        iranksum=0
        nsum=0
        irank_rep=ndim+1
        ia_rep=0
        ia_start=1
        do ia=1,ndim
          write (line,1017) ia,index2d(ia),
     -      (imemcut(index2d(ia),lev),lev=1,ncutoff)
          len=18+5*ncutoff
          line(len+1:len+1)=' '
          len=len+1
          rankav=0.0
          if (ia .eq. ndim) then
            if (nsum .gt. 0) rankav=float(iranksum)/float(nsum)
            iend=1
          else
c           Find out if cluster ends at level lev_avg
            iend=0
            do lev=1,lev_avg
              if (imemcut(index2d(ia+1),lev) .ne.
     -            imemcut(index2d(ia),lev)) iend=1
            end do
            if (index2d(ia) .lt. irank_rep) then
              irank_rep=index2d(ia)
              ia_rep=ia
            end if
            if (iend .eq. 1) then
              if (nsum .gt. 0)
     -          rankav=float(iranksum+index2d(ia))/float(nsum+1)
              iranksum=0
              nsum=0
              irank_rep=ndim+1
            else
              iranksum=iranksum+index2d(ia)
              nsum=nsum+1
            end if
          end if
          if (ilabel2d .gt. 0) then
            call lastchar(label2d(index2d(ia)),lc,llabel2d)
            line(len+1:len+lc)=label2d(index2d(ia))(1:lc)
            len=len+lc
          end if
          if (memfiles .eq. 1 .and. iend .eq. 1) then
c           Write member list file
            lmemfilename=lmemdir
            do lev=1,lev_avg
              call writeint(memfilename,lmemfilename+1,
     -          imemcut(index2d(ia),lev),ndig)
              lmemfilename=lmemfilename+ndig+1
              memfilename(lmemfilename:lmemfilename)='.'
            end do
            memfilename(lmemfilename+1:lmemfilename+4)='clst'
            lmemfilename=lmemfilename+4
c           call openfile(iu_clst,0,' ',1,'new',memfilename,
c    -        lmemfilename,notfound,0,1,1,1,0)
            open(unit=iu_clst,form='formatted',status='new',
     -        file=memfilename(1:lmemfilename),iostat=iopen)
            if (iopen .eq. 0) then
c             Write list
              do imem=ia_start,ia
                call lastchar(label2d(index2d(imem)),lc,llabel2d)
                write (iu_clst,1000) label2d(index2d(imem))(1:lc),
     -            index2d(imem)
              end do
              close (iu_clst)
            else
              write (6,*) 'ERROR: could not open ',
     -          memfilename(1:lmemfilename)
               write (6,*) 'lmemfilename=',lmemfilename
            end if
            if (iend .eq. 1) ia_start=ia+1
          end if
          if (rankav .gt. 0.0) then
            write (line(len+1:len+18),1010) rankav
            len=len+13
            call laststring(label2d(index2d(ia_rep)),ifc,ilc,lc,500)
            write (line(len+1:len+ilc-ifc+7),1022)
     -        label2d(index2d(ia_rep))(ifc:ilc)
            len=len+ilc-ifc+7
          else if (iend .eq. 1) then
            write (line(len+1:len+10),1025) index2d(ia_rep)
            len=len+10
          end if
          write (iout,1018) line(1:len)
        end do
        return
      end if
      isubcl=0
      if (nclust .lt. ndim .and. iclstyp .eq. 1) then
c       Sub-clustering only works for single-link clustering
        call askyn('Do you want to try sub clustering',33,1,-1,isubcl,
     -    29,0)
        if (isubcl .gt. 0) then
          call quiz(ans,isubclustertyp,' ',' ',0,
     -      'subclustering algorithm',23,0,5,6,78)
          call getint('Minimum number of members for subclustering',43,
     -      999999,1,ndim,minsubclust,30)
          write (iout,1011) minsubclust
          if (isubclustertyp .eq. 1) then
            write (iout,1012)
          else if (isubclustertyp .eq. 2) then
            call getreal(
     -        'Maximum percent of density for being on the surface',51,
     -        20.0,surfacepercent,1,79)
            surfacefract=surfacepercent/100.0
            call getreal(
     -        'Distance exponent for local density descriptor',46,
     -        2.0,denexp,0,81)
            denexphalf=denexp/2.0
            write (iout,1013) surfacepercent,denexp
          end if
          do ic=nclust,1,-1
            if (ilclst(ic)-ifclst(ic) .gt. minsubclust) then
              call rmsdsubcluster(ic,ifclst(ic),ilclst(ic),index2d,
     -          ixclst,ifa_s,ila_s,rdclust,isubclustertyp,
     -          surfacefract,denexphalf,it1,it2,indexa,ih,it4,value,
     -          cv,iout,nclustic,label,llabel)
              if (nclustic .gt. 1) then
                write (iout,2065) ic,nclustic
                write (6,2065) ic,nclustic
c               Adjust cluster limits
                do icc=nclust,ic+1,-1
                  ifclst(icc)=ifclst(icc-(nclustic-1))
                  ilclst(icc)=ilclst(icc-(nclustic-1))
                end do
                do icc=1,nclustic
                  ifclst(ic+icc-1)=ifa_s(icc)
                  ilclst(ic+icc-1)=ila_s(icc)
                end do
                call sortlist(iout,index2d,ilclst(nclust),it1,it2,'IX2',
     -            1,mx2d)
                call reportclust(ndim,ic,ic,ic+nclustic-1,ifclst,ilclst,
     -            index2d,value,it1,ifa_s,ila_s,ih,cv,indexa,irepav,
     -            irepmx,irepeng,irepkm,engcl,nhbdist,etotsaved,
     -            ietotsaved,ifindbestrep,label,llabel,isorttype,
     -            idistprint,nomemprint,iout,maxframe,mx2d)
                print *,'ic=',ic,' IREPMX=',irepmx(ic)
              else
                write  (iout,*) 'None of the clusters split up'
                write  (6,*) 'None of the clusters split up'
              end if
            end if
          end do
        end if
      end if
      if (ilabel2d .eq. 1 .and. iclstyp .ne. 3) then
        write (iout,1021)
        do ic=1,nrdclust
          do i=ifclst(ic),ilclst(ic)
            lmemlab=0
            if (i .eq. ifclst(ic)) then
              memlab(1:6)=' FIRST'
              lmemlab=6
            end if
            if (ixclst(i) .eq. irepav(ic)) then
              memlab(lmemlab+1:lmemlab+6)=' AVMIN'
              lmemlab=lmemlab+6
            end if
            if (ixclst(i) .eq. irepmx(ic)) then
              memlab(lmemlab+1:lmemlab+6)=' MXMIN'
              lmemlab=lmemlab+6
            end if
            call lastchar(label2d(ixclst(i)),llab,llabel2d)
            write (iout,1020) ic,ixclst(i),label2d(ixclst(i))(1:llab),
     -        memlab(1:lmemlab)
          end do
        end do
      end if
      if (ietotsaved .gt. 0) then
        call askyn(
     -    'Do you want a histogram of cluster member energies',50,1,-1,
     -    ihiste,0,0)
        if (ihiste .gt. 0)
     -    call getreal('Histogram bin size',18,1.0,rbin,1,000)
        call askyn(
     -    'Do you want to list the energies of the cluster members',55,
     -    1,-1,iliste,0,0)
        if (ihiste+iliste .gt. 0) then
          do ic=1,nclust
            do i=ifclst(ic),ilclst(ic)
              cv(i-ifclst(ic)+1)=etotsaved(1,index2d(i))
            end do
            nmem=ilclst(ic)-ifclst(ic)+1
            if (iliste .gt. 0) write (iout,1026) ic,(cv(i),i=1,nmem)
            if (ihiste .gt. 0) call histogram(cv,nmem,0.0,rbin,
     -        'cluster energies',16,it1,iout,mx2d)
          end do
        end if
      end if
      call askyn('Do you want to run more clustering',34,1,-1,
     -  ians,0,0)
      if (ians .eq. 1) then
        call trnsfi(index2d,it2d,ndim)
        go to 8910
      end if
      return
1000  format(a,i8)
1010  format(' <R>=',f8.1)
1011  format(' Subclustering clusters containing at least',i4,
     -  ' members')
1012  format(' Subclustering ignores members with the ',
     -  'smallest number of neighbors')
1013  format(' Subclustering ignores members whose local density is ',
     -  'less than ',f5.1,' %',/,' of the maximum local density',/,
     -  ' Local density is the sum of 1/rij^',f3.1)
1014  format(a,' cutoff for clustering')
1015  format('Cutoff # ',i2)
1016  format(' Cutoff list: ',5f10.5)
1017  format(i6,' (',i6,') : ',5(i4,'.'))
1018  format(1x,a,1x,a)
1019  format(/,' In the list below, each line corresponds to one ',
     -  'clustered item',/,' For  ( N ) i1. i2. i3. i4. i5., N is the ',
     -  ' original index of this item,'/,' i1 is the cluster number N ',
     -  ' belongs to, obtained with cutoff #1; ',/,' i2 is the ',
     -  'cluster number of N within cluster i1, obtained with cutoff ',
     -  '#2;',/,' i3 is the cluster number of N within cluster i2, ',
     -  'obtained with cutoff #3; etc.',/
     -  ' <R> is the average rank of the cluster members',/)
1020  format(' Cluster #',i5,' member #',i5,1x,a,a)
1021  format(/,' List of cluster member labels read from the distance ',
     -  'matrix file')
1022  format(' REP: ',a)
1023  format(' Cluster rank averages and representatives are obtained ',
     -  'at cutoff level',i2)
1024  format(' Directory ',a,' is not found',/,
     -  ' Create the directory and try again')
1025  format(' S R=',i5)
1026  format(' Cluster #',i6,' Energies:',/,(5e12.5))
2060  format(/,' === Clustering method: ',a,/,
     -  ' Number of items to cluster=',i6)
2061  format(' Number of clusters requested=',i5)
2062  format(' Clusters are reported in the order of increasing ',
     -  'average indices (<rank>)')
2065  format(' Cluster ',i4,' partitioned into',i3,' subclusters')
2066  format(1x,a,' threshold resulting in ',i3,' clusters=',f6.2)
2090  format(1x,a,' threshold=',f9.2,' number of clusters=',i4)
2092  format(' Clustering with ',a,' threshold=',f8.3,' A',/)
2093  format(' Clustering into ',i4,' sets. Range of RMSDs: [',
     -  f8.2,',',f8.2,']'/)
2094  format(' Clustering failure: instead of',i4,' clusters',i4,
     -  ' was generated')
      end
