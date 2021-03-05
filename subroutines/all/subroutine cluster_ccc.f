      subroutine cluster_ccc(ndials,angnames,langnames,
     -  ixclst,ifclst,ilclst,index2d,it1,it2,it3,it4,label,llabel,iout,
     -  inpfile,namleni)
      dimension langnames(ndials),ixclst(ndials),ifclst(ndials),
     -  ilclst(ndials),index2d(ndials),it1(ndials),it2(ndials),
     -  it3(ndials),it4(ndials)
      character*(*) angnames(ndials),label
      character*(*) inpfile
      parameter (MAXPHI=400,MAX2D=5000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-(MAX2D*MAX2D+9*MAX2D))
      common /nnwork/ ccc(MAX2D,MAX2D),irepav(MAX2D),irepmx(MAX2D),
     -  irepeng(MAX2D),irepkm(MAX2D),etotsaved(2,MAX2D),cv(MAX2D),
     -  engcl(MAX2D),indexa(MAX2D),fill(IFILL2)
      character*41 clstyp
      common /cluster_typ/ nclstyp,inumclst(9),ireadcutoff(9),
     -  lclstyp(9),clstyp(9)
      character*1 ans
      character*80 line
      do i=1,ndials
        do j=1,ndials
          ccc(i,j)=1.0-abs(ccc(i,j))
        end do
      end do
      nrdclust=0
8910  call quiz(ans,iclstyp,' ',' ',0,'clustering algorithm',20,
     -  0,5,6,40)
      if (iclstyp .gt. 7) then
        print *,'Sorry, this option is not available here'
        go to 8910
      end if
      write (iout,2060) clstyp(iclstyp),ndials
      write (6,2060) clstyp(iclstyp),ndials
      if (iclstyp .ne. 4) call zeroiti(irepkm,0,ndials)
      if (inumclst(iclstyp) .eq. 1) then
        call getint('Number of clusters requested',28,999999,1,ndials,
     -    nclust,0)
        rdclust=0.0
        write (iout,2061) nclust
      else if (ireadcutoff(iclstyp) .eq. 1) then
        write (line,1014) label(1:llabel)
        call getreal(line,22+llabel,999999.0,rdclust,1,0)
      end if
      call rmsdcluster(rdclust,1,ndials,index2d,iwt,ixclst,ifclst,
     -    ilclst,it1,it2,it3,irepav,irepmx,it4,c,cent,cent_prev,irepkm,
     -    0,iclstyp,0,0,nclust,1,0,ifail,1,' ',1,iout)
      write (iout,7824) (ifclst(i),ilclst(i),i=1,nclust)
 7824 format(' After RMSDCLUSTER Cluster limits: ',('[',i5,',',i5,']'))
      write (iout,6792) 'IXCLST',(ixclst(i),i=1,ndials)
 6792 format(' AFTER RMSDCLUSTER ',a,':',/,(20i5))
      if (ifail .gt. 0) go to 8910
c     Members of cluster ic: (index2d(i),i=ifclst(ic),ilclst(ic))
      if (ans .ne. 'k') write (6,2090) label(1:llabel),rdclust,nclust
      call reportclust(ndials,0,1,nclust,ifclst,ilclst,index2d,value,
     -  it1,it2,it3,it4,cv,indexa,irepav,irepmx,irepeng,irepkm,engcl,
     -  nhbdist,etotsaved,0,1,'RMSD',4,1,idistprint,nomemprint,iout,
     -  MAX2D,MAX2D)
      return
1014  format(a,' cutoff for clustering')
2060  format(/,' === Clustering method: ',a,/,
     -  ' Number of items to cluster=',i6)
2061  format(' Number of clusters requested=',i5)
2090  format(1x,a,' threshold=',f9.2,' number of clusters=',i4)
      end
