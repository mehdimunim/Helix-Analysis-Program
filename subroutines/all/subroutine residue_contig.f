      subroutine residue_contig(n,ires,isegno,index,ixseg,isr,ifa,ila,
     -  itmp1,itmp2,maxres,max)
      dimension ires(max),isegno(max),index(max),ixseg(maxres),isr(max),
     -  ifa(max),ila(max),itmp1(max),itmp2(max)
c     print *,'RESIDUE_CONTIG n=',n
      call indexit(index,1,n,0)
      do i=1,n
        isr(i)=(max+1)*ixseg(isegno(i))+ires(i)+1
      end do
      call mrgsrti(6,index,isr,n,ifa,ila,itmp1,itmp2,max)
      return
      end
