      subroutine printclusterpdb(c,frocc,bfac,ncl,ifclst,ilclst,ixclst,
     -  title,na)
      dimension c(3,na),frocc(na),bfac(na),ifclst(ncl),ilclst(ncl),
     -  ixclst(na)
      character*80 title
      character*200 pdbout
c     Print a PDB file of atoms where each cluster is a separate residue
c     print *,'PRINTCLUSTERPDB ncl=',ncl,' na=',na
      nl=0
      call openfile(50,0,'clustered PDB',13,'new',pdbout,nl,notfnd,
     -  0,1,1,0,0)
      write (50,1001) 'REMARK '//title(1:72)
      write (50,1001) 'REMARK Atoms clustered into different residues'
      na=1
      write (50,1000) na,' O  ','CNT','A',1,0.0,0.0,0.0,1.0,0.0
      do icl=1,ncl
        do ia=ifclst(icl),ilclst(icl)
          na=na+1
          iap=ixclst(ia)
          write (50,1000) na,' H  ','CLS','C',icl,
     -      (c(k,iap),k=1,3),frocc(iap),bfac(iap)
        end do
      end do
      write (50,1001) 'END'
      close (50)
      return
1000  format('ATOM  ',i5,1x,a4,1x,a3,1x,a1,i4,1x,3x,3f8.3,2f6.2)
1001  format(a)
      end
