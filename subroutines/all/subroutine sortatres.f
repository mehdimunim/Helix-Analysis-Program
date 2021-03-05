      subroutine sortatres(line,n,index,indexn,indexo,ifa,ila,iresn,
     -  ireso,isegn,iresncol1,iresncol2,isegno,marker,inpcrdtyp,nconfig,
     -  iout,maxrepconf,maxrec)
      character* 132 line(maxrec)
      character*6 marker
      dimension index(n),isegno(n),indexo(n),indexn(n),ifa(n),ila(n),
     -  iresn(n),ireso(n),isegn(n)
      character*5 crdext
      common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
     -  iommod,iommc,iommc4,iogro,iomol2,iomae,iocif,ioxxx,ioins,
     -  ionxyz,iosxyz,iosxyzrq,iograsp,iofull,lext(19),crdext(19)
c     Rearrange coordinates to keep atoms of the same residue together
c     print *,'SORTATRES,iresncol1,2',n,iresncol1,iresncol2
      if (n .eq. 1) return
      if (nconfig .le. maxrepconf)
     -  print *,'Sorting by segment (chain) id and residue id'
      if (n .gt. 1000) print *,'..... Wait'
      nresncol=iresncol2-iresncol1+1
      call indexit(indexo,1,n,0)
c     Individual intervals consist of single elements at start
      nn=n
      call indexit(ifa,1,n,0)
      call indexit(ila,1,n,0)
c     Merge pairs of intervals
11    nnpair=nn/2
      l=1
      do i=1,nnpair
        call mergerec(line,index,indexo,isegno,iresncol1,iresncol2,
     -    ifa(l),ila(l),ifa(l+1),ila(l+1),iresn,ireso,isegn,n,maxrec)
        ifa(i)=ifa(l)
        ila(i)=ila(l+1)
        l=l+2
      end do
      if (2*nnpair .eq. nn) go to 21
c     Take care of last (odd) interval
      call mergerec(line,index,indexo,isegno,iresncol1,iresncol2,
     -  ifa(nnpair),ila(nnpair),ifa(nn),ila(nn),iresn,ireso,isegn,n,
     -  maxrec)
      ila(nnpair)=ila(nn)
21    nn=nnpair
      if (nn .gt. 1) go to 11
      do i=1,n
        indexn(indexo(i))=i
      end do
      if (inpcrdtyp .eq. iommod)
     -  call nnupdate(index,indexn,indexo,n,n,line,0,maxrec)
      if (marker .ne. '      ') write (iout,1270) marker
      return
1270  format(a6,' Atoms sorted by residues by SIMULAID')
      end
