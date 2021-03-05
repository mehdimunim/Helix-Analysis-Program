      subroutine atomdel(line,idcol,asterisk,n,nslt,nrecdel,c,isegno,
     -  index,indexn,indexo,inpcrdtyp,iatnum,charge,altcol,inscol,
     -  iresno,ixres,atnames,ifres,ilres,nconfig,ndel,maxrsd,maxrec)
      character*132 line(maxrec)
      character*1 altcol(maxrec),inscol(maxrec)
      character*8 atnames(maxrec)
      dimension c(3,n),isegno(n),index(n),indexn(n),indexo(n),iatnum(n),
     -  iresno(n),ixres(n),charge(n),ifres(maxrsd),ilres(maxrsd)
      character*1 asterisk
      character*5 crdext
      common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
     -  iommod,iommc,iommc4,iogro,iomol2,iomae,iocif,ioxxx,ioins,
     -  ionxyz,iosxyz,iosxyzrq,iograsp,iofull,lext(19),crdext(19)
c     print *,'ATOMDEL nrecdel,n,idcol,nconfig=',nrecdel,n,idcol,nconfig
      if (nrecdel .eq. 0) return
      ndel=0
      nsltdel=0
      do i=1,n
c       write (77,*) 'ATOMDEL i=',i,' AC:',line(index(i))(idcol:idcol)
c       Update index to eliminate the dropped molecules
        if (line(index(i))(idcol:idcol) .eq. asterisk) then
          ndel=ndel+1
          if (i .le. nslt) nsltdel=nsltdel+1
          indexo(i)=0
        else if (ndel .gt. 0) then
          index(i-ndel)=index(i)
          indexo(i)=i-ndel
          iatnum(i-ndel)=iatnum(i)
          iresno(i-ndel)=iresno(i)
          isegno(i-ndel)=isegno(i)
          charge(i-ndel)=charge(i)
          call trnsfr(c(1,i-ndel),c(1,i),3)
          altcol(i-ndel)=altcol(i)
          inscol(i-ndel)=inscol(i)
          atnames(i-ndel)=atnames(i)
        end if
      end do
      if (nsltdel .gt. 0 .and. nconfig .eq. 0)
     -  print *,'Number of solute atoms deleted=',nsltdel
      if (ndel .ne. nrecdel) then
        write (6,1000) asterisk,ndel,nrecdel
      end if
c     print *,'atomdel n,inpcrdtyp=',n,inpcrdtyp
      if (inpcrdtyp .eq. iommod)
     -  call nnupdate(index,indexn,indexo,n-ndel,n,line,1,maxrec)
      n=n-ndel
      nslt=nslt-nsltdel
c     Establish atom ranges for the residues reflecting the deletions
      ifres(1)=1
      numres=0
      ixres(1)=1
      do ia=2,n
        if (iresno(ia) .ne. iresno(ia-1) .or.
     -      isegno(ia) .ne. isegno(ia-1)) then
          numres=numres+1
          ilres(numres)=ia-1
          do ja=ifres(numres),ilres(numres)
            ixres(ja)=numres
          end do
          ifres(numres+1)=ia
        end if
      end do
      ilres(numres+1)=n
      nrecdel=0
      return
1000  format(' PROGRAM ERROR: Number of ',a,'-marked lines found=',
     -  i9,/,5x,'differs from number of records to be deleted=',i9)
      end
