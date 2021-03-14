      subroutine comparetree(nneig,ineig,ian,ia,in1,in2,idiff,n,maxneig)
      dimension ian(n),nneig(n),ineig(maxneig,n)
      dimension iused1(1000),iused2(1000),nnpar1(1000),nnoffsp1(1000),
     -  nnpar2(1000),nnoffsp2(1000),nian1(100),nian2(100)
      data noffsp1 /0/,noffsp2 /0/
c     print *,'COMPARETREE ia,in1,in2=',ia,in1,in2
      call zeroiti(nian1,0,100)
      call zeroiti(nian2,0,100)
      call zeroiti(iused1,0,n)
      call zeroiti(iused2,0,n)
      iused1(ia)=1
      iused1(in1)=1
      iused2(ia)=1
      iused1(in2)=1
      npar1=1
      nnpar1(1)=in1
      npar2=1
      nnpar2(1)=in2
      lev=0
      idiff=0
      do while ((lev .eq. 0 .or. max0(noffsp1,noffsp2) .gt. 0) .and.
     -          idiff .eq. 0)
        noffsp1=0
        do ip=1,npar1
          do in=1,nneig(nnpar1(ip))
            ja=ineig(in,nnpar1(ip))
            if (iused1(ja) .eq. 0) then
              noffsp1=noffsp1+1
              nnoffsp1(noffsp1)=ja
              nian1(ian(ja))=nian1(ian(ja))+1
              iused1(ja)=1
            end if
          end do
        end do
        noffsp2=0
        do ip=1,npar2
          do in=1,nneig(nnpar2(ip))
            ja=ineig(in,nnpar2(ip))
            if (iused2(ja) .eq. 0) then
              noffsp2=noffsp2+1
              nnoffsp2(noffsp2)=ja
              nian2(ian(ja))=nian2(ian(ja))+1
              iused2(ja)=1
            end if
          end do
        end do
        idiff=0
        do i=1,50
          idiff=idiff+iabs(nian2(i)-nian1(i))
        end do
c        write (6,7677) 'offsp1',(nnoffsp1(i),i=1,noffsp1)
c        write (6,7677) 'offsp2',(nnoffsp2(i),i=1,noffsp2)
c7677    format(1x,a,':',(30i3))
        if (idiff .eq. 0) then
          npar1=noffsp1
          if (npar1 .gt. 0) call trnsfi(nnpar1,nnoffsp1,npar1)
          npar2=noffsp2
          if (npar2 .gt. 0) call trnsfi(nnpar2,nnoffsp2,npar2)
        end if
        lev=lev+1
      end do
      return
      end
