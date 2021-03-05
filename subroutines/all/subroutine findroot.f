      subroutine findroot(lev,ia,ianoneig,iatno,nneig,ineig,line,index,
     -  inpcrdtyp,iofull,ia1,n,ifail,maxneig,maxrec)
      dimension index(n),iatno(n),nneig(n),ineig(maxneig,n)
      character*2 lev
      character* 132 line(maxrec)
      character*16 lab1,lab2
      dimension iahv(8)
      data lab1 /'Possible R2 atom'/,lab2 /'Index of atom R2'/,
     -  ia2heavy /0/
      if (ifail .gt. 0) return
      lab1(10:11)=lev
      lab2(15:16)=lev
c     print *,'FINDROOT ia=',ia,' inpcrdtyp=',inpcrdtyp,' ix(ia)=',index(ia)
      nnheavy=0
      if (nneig(ia) .eq. 0) then
        write (6,1001) ia,' '
        ifail=1
        return
      end if
      do ja=1,nneig(ia)
        ia1=ineig(ja,ia)
c       print *,'ia,ja,ia1,iatno(ia1)=',ia,ja,ia1,iatno(ia1)
        if (iatno(ia1) .gt. 1 .and. ia1 .ne. ianoneig) then
          call listatom(line,index,iatno,ia1,inpcrdtyp,iofull,lab1,16,n,
     -      maxrec)
          nnheavy=nnheavy+1
          iahv(nnheavy)=ia1
          ia2heavy=ia1
        end if
      end do
      if (nnheavy .gt. 1) then
100     call getint(lab2,16,0,1,n,ia1,46)
        nfound=0
        do ja=1,nnheavy
          if (ia1 .eq. iahv(ja)) nfound=1
        end do
        if (nfound .eq. 0) then
          write (6,1000) ia1,(iahv(ja),ja=1,nnheavy)
          go to 100
        end if
      else if (nnheavy .eq. 1) then
        ia1=ia2heavy
        write (6,1002) ia1,lev
      else
        write (6,1001) ia,' heavy atom '
        ifail=1
        return
      end if
      return
1000  format(' ERROR: index',i6,' is not among indices ',8i6)
1001  format(' ERROR: atom',i6,' has no',a,'neighbor')
1002  format(' Atom',i6,' will be used for ',a)
      end
