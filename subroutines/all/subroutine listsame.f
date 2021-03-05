      subroutine listsame(ineig,ian,nn,n,nnpairs,nsamepair)
      dimension ineig(nn),ian(n),nnpairs(2,6)
      nsamepair=0
      do ia=2,4
        do ja=1,ia-1
          if (ian(ineig(ia)) .eq. ian(ineig(ja))) then
            nsamepair=nsamepair+1
            nnpairs(1,nsamepair)=ineig(ia)
            nnpairs(2,nsamepair)=ineig(ja)
          end if
        end do
      end do
      return
      end
