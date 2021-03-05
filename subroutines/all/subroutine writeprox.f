      subroutine writeprox(iout,irarep,inarep,line,index,r2,iresno,
     -  ir1,ir2,ic1,ic2,is1,is2,maxrec)
      dimension index(maxrec),iresno(maxrec)
      character* 132 line(maxrec)
c      write (6,1000) index(irarep),index(inarep)
c1000  format(' index(ir/n)=',2i4)
      r2=sqrt(r2)
      write (iout,2060) line(index(irarep))(is1:is2),
     -  line(index(irarep))(ir1:ir2),iresno(irarep),
     -  line(index(irarep))(ic1:ic2),
     -  irarep,line(index(inarep))(is1:is2),
     -  line(index(inarep))(ir1:ir2),iresno(inarep),
     -  line(index(inarep))(ic1:ic2),inarep,r2
      return
2060  format(1x,a,1x,a,' (',i4,') ',a,' (',i5,') - ',a,1x,a,
     -  ' (',i4,') ',a,' (',i5,')  Dist=',f8.3,' A')
      end
