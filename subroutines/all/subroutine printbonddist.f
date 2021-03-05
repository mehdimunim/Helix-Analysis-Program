      subroutine printbonddist(nframe,nbfound,nhbdist,index2d,
     -  ianc_anc,naabond,iout,bondtype,lbondtype,iselfanc,maxbondcount,
     -  maxhb)
      dimension nhbdist(maxhb),index2d(maxhb),
     -  ianc_anc(maxhb)
      character*(*) bondtype
      dimension nhbhist(10)
      maxbondcount=0
      minbondcount=nframe
      nbondsum=0
      naabond=0
      do ib=1,nbfound
        nhbd=nhbdist(index2d(ib))
        if (nhbd .gt. maxbondcount) maxbondcount=nhbd
        if (nhbd .lt. minbondcount) minbondcount=nhbd
        nbondsum=nbondsum+nhbd
        naabond=naabond+ianc_anc(ib)
      end do
      call zeroiti(nhbhist,0,10)
      do ib=1,nbfound
        ix=(10*nhbdist(index2d(ib))-1)/maxbondcount+1
        nhbhist(ix)=nhbhist(ix)+1
      end do
      write (6,2000) bondtype(1:lbondtype),nbfound,maxbondcount
      write (6,2001) bondtype(1:lbondtype),naabond
      if (iselfanc .eq. 1)
     -  write (iout,2001) bondtype(1:lbondtype),naabond
      write (6,2003) bondtype(1:lbondtype),float(nbondsum)/float(nframe)
      write (iout,2003) bondtype(1:lbondtype),
     -  float(nbondsum)/float(nframe)
      write (6,2002) nframe,minbondcount,maxbondcount,
     -  (nhbhist(i),i=1,10)
      write (iout,2002) nframe,minbondcount,maxbondcount,
     -  (nhbhist(i),i=1,10)
      return
2000  format(' Number of ',a,'-bond types found=',i5,/,
     -  ' Most persistent bond was present in ',i6,' configurations')
2001  format(' Number of anchor-anchor ',a,' bonds found=',i6)
2002  format(' Histogram of the number of frames a bond is present ',
     -  'over',i7,' configurations',/,
     -  ' (10% intervals of [',i5,' - ',i5,'] range:',/,10i7)
2003  format(' Total number of ',a,' bonds / number of frames=',f7.1)
      end
