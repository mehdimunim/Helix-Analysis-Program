      subroutine maybehbond(r2,i,j,nneig,nhbneig,ineig,inamcol1,
     -  inamcol2,rlimhb2,atnami,resnami,n,maxng)
      dimension nneig(n),nhbneig(n),ineig(maxng,n)
      character*4 atnami
      character*8 resnami
c     write (78,*) 'MAYBEHBOND i,j=',i,j
      if (r2 .lt. rlimhb2) then
        maxngij=max0(nneig(i)+nhbneig(i),nneig(j)+nhbneig(j))
        if (maxngij .lt. maxng) then
c         Hydrogen bond found
          ineig(maxng-nhbneig(i),i)=j
          nhbneig(i)=nhbneig(i)+1
          ineig(maxng-nhbneig(j),j)=i
          nhbneig(j)=nhbneig(j)+1
        else
          if (inamcol2 .lt. inamcol1) then
            write (6,2000) i,maxng
          else
             write (6,2001) i,atnami,resnami,maxng
          end if
          write (6,2002) maxngij
        end if
      end if
      return
2000  format(' ERROR: atom',i6,' has more than ',i2,' neighbours',
     -  ' and hydrogen bonds')
2001  format(' ERROR: atom',i6,' (',a4,a6,') has more than ',i2,
     -  ' neighbours and hydrogen bonds')
2002  format(' - increase the parameter MAXNEIG and recompile')
      end
