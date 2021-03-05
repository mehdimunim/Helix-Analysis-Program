      subroutine print1charlist(charlist,nlist,line,index,ifres,
     -  irescol1,irescol2,maxrec)
      character*1 charlist(nlist)
      character*(*) line(maxrec)
      dimension index(maxrec),ifres(maxrec)
      character*1 aa1(50)
      character*8 resname
      lresname=irescol2-irescol1+1
      i0=0
      do while (i0 .lt. nlist)
        print *
        imax=i0+min0(50,nlist-i0)
        write (6,1001)'Chirality:',(charlist(i),i=i0+1,imax)
c       Get 1-character AA list
        iimax=imax-i0
        do i=i0+1,imax
          ii=i-i0
          resname(1:lresname)=line(index(ifres(i)))(irescol1:irescol2)
          call leftadjustline(resname,1,lresname)
          call changeprot(resname,aa1(ii),2)
c          if (i0 .eq. 0)
c     -      write (6,7711) i,ii,index(i),lresname,resname,aa1(ii)
c7711  format('  i,ii=',2i4,' index=',i4,' lresname=',i4,
c     -  ' resname=',a,' aa1=',a)
        end do
        write (6,1001)'Residue:  ',(aa1(ii),ii=1,iimax)
        moddiv=1
        nline=alog10(float(imax))+1
        do iline=1,nline
          write (6,1002) (mod(i/moddiv,10),i=i0+1,imax)
          moddiv=moddiv*10
        end do
        i0=i0+50
      end do
      return
1001  format(1x,a10,50a1)
1002  format(11x,50i1)
      end
