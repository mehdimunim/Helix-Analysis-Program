      subroutine hph_sltblist(ibondtype,c,n,nhbneig,ineig,iout,ic1,ic2,
     -  in1,in2,ir1,ir2,irn1,irn2,line,index,bondname,lbondname,rlab,
     -  lrlab,maxng)
      dimension c(3,n),ineig(maxng,n),nhbneig(n),index(n),lbondname(5)
      character*132 line(n)
      character*22 bondname(5)
      character*(*) rlab
c     print *,'HPHLIST n,maxng=',n,maxng
      write (iout,1000) bondname(ibondtype)(1:lbondname(ibondtype))
      do ia=1,n
        do jaa=1,nhbneig(ia)
          ja=ineig(maxng-jaa+1,ia)
          call readint(line(index(ia)),irn1,irn2,irnia,2,1,irerr)
          call readint(line(index(ja)),irn1,irn2,irnja,2,1,irerr)
          write (iout,1001) line(index(ia))(ic1:ic2),
     -      ia,line(index(ia))(in1:in2),irnia,line(index(ia))(ir1:ir2),
     -      line(index(ja))(ic1:ic2),ja,line(index(ja))(in1:in2),
     -      irnja,line(index(ja))(ir1:ir2),rlab(1:lrlab),
     -      sqrt(dist2(c(1,ia),c(1,ja)))
        end do
      end do
      return
1000  format(' List of ',a,'(s) found:')
1001  format(1x,a,i5,1x,a4,i4,1x,a4,' - ',a,1x,i5,1x,a4,i4,1x,a4,
     -  ' r(',a,')=',f6.1)
      end
