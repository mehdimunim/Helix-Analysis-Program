      subroutine nnlist0(nfirst,n,nslt,islvw,iatnum,ifchrg,c,nneig,
     -  nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,npneig,line,
     -  irescol1,irescol2,inamcol1,inamcol2,index,maxng,hblimfac,angmin,
     -  ihbondcalc,indices,nbox,ixres,isegno,ifail,nframe,radtodeg,
     -  maxbox,maxrec,LEVTEST)
      dimension nneig(n),ineig(maxng,n),iatnum(n),ifchrg(n),c(3,n),
     -  nhbneig(n),nhneig(n),nnneig(n),ncneig(n),nsneig(n),npneig(n),
     -  index(n),indices(maxbox,maxrec),nbox(maxrec),ixres(n),isegno(n)
      common /connatdat/ ramax(99),ramax2(99),hlimfac,ianfg(99),
     -  namfcg(100),nrmw
      character* 132 line(maxrec)
c     print *,'NNL0 n,nslt,islvw,ihbondcalc=',n,nslt,islvw,ihbondcalc
      ifail=0
      if (n .lt. nfirst) return
      if (n .eq. 1) then
        nneig(1)=0
        nhbneig(1)=0
        return
      end if
      nn=0
      do ia=nfirst,n
        if (isegno(ia) .eq. -1) nn=nn+1
      end do
      if (nn .gt. 0 .or. ihbondcalc .eq. 1) then
c       LES structure, call nnlist00 with the whole range
        call nnlist00(nfirst,n,nslt,islvw,iatnum,ifchrg,c,nneig,nhbneig,
     -    ineig,nhneig,nnneig,ncneig,nsneig,npneig,line,irescol1,
     -    irescol2,inamcol1,inamcol2,index,maxng,hblimfac,angmin,
     -    indices,nbox,ixres,isegno,ifail,nframe,radtodeg,maxbox,maxrec,
     -    LEVTEST)
        if (ifail .gt. 0) go to 2000
        return
      else
c       Call the actual near neighbour search by segments
        isg=isegno(nfirst)
        ifat=nfirst
        ilat=nfirst
        do while (isg .le. isegno(n))
c         Find limits of the next segment
          do while (isegno(ilat) .eq. isg .and. ilat .lt. n)
            ilat=ilat+1
          end do
          if (isegno(ilat-1) .ne. isegno(ilat) .and. ilat .le. n)
     -       ilat=ilat-1
          if (ilat .lt. n) then
c           If segments are not contiguous, just call nnlist00 for
c           the whole range
            if (isegno(ilat+1) .ne. isegno(ilat)+1) then
              write (6,1000) ilat,isegno(ilat),ilat+1,isegno(ilat+1)
              call nninit(nneig,nhbneig,ineig,nhneig,nnneig,ncneig,
     -          nsneig,npneig,nfirst,n,1,maxng)
              call nnlist00(nfirst,n,nslt,islvw,iatnum,ifchrg,c,nneig,
     -          nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,npneig,line,
     -          irescol1,irescol2,inamcol1,inamcol2,index,maxng,
     -          hblimfac,angmin,indices,nbox,ixres,isegno,ifail,nframe,
     -      radtodeg,maxbox,maxrec,LEVTEST)
              if (ifail .gt. 0) go to 2000
              return
            end if
          end if
          isg=isg+1
          call nnlist00(ifat,ilat,nslt,islvw,iatnum,ifchrg,c,nneig,
     -      nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,npneig,line,
     -      irescol1,irescol2,inamcol1,inamcol2,index,maxng,hblimfac,
     -      angmin,indices,nbox,ixres,isegno,ifail,nframe,radtodeg,
     -      maxbox,maxrec,LEVTEST)
          ifat=ilat+1
          ilat=ilat+1
        end do
      end if
      if (ifail .eq. 0) return
2000  call nninit(nneig,nhbneig,ineig,nhneig,nnneig,ncneig,
     -  nsneig,npneig,nfirst,n,1,maxng)
      call nnlist0o(nfirst,n,iatnum,c,nneig,nhbneig,ineig,nhneig,nnneig,
     -  ncneig,nsneig,npneig,line,irescol1,irescol2,inamcol1,inamcol2,
     -  index,maxng,hblimfac,maxrec)
      return
1000  format(' Segments are not contiguous: isegno(',i5,')=',i4,
     -  ' isegno(',i5,')=',i4)
      end
