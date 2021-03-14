      subroutine nnlist0o(nfirst,n,iatnum,c,nneig,nhbneig,ineig,nhneig,
     -  nnneig,ncneig,nsneig,npneig,line,irescol1,irescol2,inamcol1,
     -  inamcol2,index,maxng,hblimfac,maxrec)
      dimension nneig(n),ineig(maxng,n),iatnum(n),c(3,n),nhbneig(n),
     -  nhneig(n),nnneig(n),ncneig(n),nsneig(n),npneig(n),index(n)
      common /connatdat/ ramax(99),ramax2(99),hlimfac,ianfg(99),
     -  namfcg(100),nrmw
      character* 132 line(maxrec)
      character*4 atnami
      character*8 resnami
c     Set up neighbour list using the trivial algorithm
c     nhneig, ncneig,nnneig,nsneig,npneig:
c     Number of H, C, N, S, P neighbours, resp.
      if (n .eq. nfirst) return
      do i=nfirst,n
        if (inamcol2 .ge. inamcol1) then
          resnami='     '
          resnami=line(index(i))(irescol1:irescol2)
          atnami=line(index(i))(inamcol1:inamcol2)
        end if
        j1=i+1
c       write (77,7711) i,iatnum(i),ramax(iatnum(i)),ramax2(iatnum(i))
c7711    format(i5,' iano=',i2,' ramax,2=',2f10.5)
        do j=j1,n
          r2=dist2(c(1,i),c(1,j))
          call decidebondcut(iatnum(i),iatnum(j),rlim)
          if (r2 .le. rlim .and. iatnum(i)+iatnum(j) .gt. 2 .and.
     -        ramax2(iatnum(i))*ramax2(iatnum(j)) .gt. 0.0) then
c           Bond found. Bonds to electron, charge or lone pair are not
c           saved in the 'real' atoms' list here.
            call savebond(j,i,iatnum,nneig,ineig,nhbneig,nhneig,ncneig,
     -        nnneig,npneig,nsneig,resnami,atnami,maxng,n,1,
     -        inamcol1,inamcol2,0)
            call savebond(i,j,iatnum,nneig,ineig,nhbneig,nhneig,ncneig,
     -        nnneig,npneig,nsneig,resnami,atnami,maxng,n,1,
     -        inamcol1,inamcol2,0)
          else if (iatnum(i)  .eq. 1 .and. iatnum(j) .gt. 1 .and.
     -             iatnum(j) .ne. 6 .or.
     -             iatnum(j)  .eq. 1 .and. iatnum(i) .gt. 1 .and.
     -             iatnum(i) .ne. 6 )  then
c           No bond, check if hydrogen bond
            if (r2 .lt. rlim*(hblimfac/hlimfac)) then
              if (max0(nneig(i)+nhbneig(i),nneig(j)+nhbneig(j)) .lt.
     -            maxng) then
c               Hydrogen bond found
                ineig(maxng-nhbneig(i),i)=j
                nhbneig(i)=nhbneig(i)+1
                ineig(maxng-nhbneig(j),j)=i
                nhbneig(j)=nhbneig(j)+1
              else
                if (inamcol2 .lt. inamcol1)
     -            write (6,2000) i,maxng
                if (inamcol2 .ge. inamcol1)
     -            write (6,2001) i,atnami,resnami,maxng
                write (6,2002)
              end if
            end if
          end if
        end do
      end do
      return
2000  format(' ERROR: atom',i6,' has more than ',i2,' neighbours',
     -  ' and hydrogen bonds')
2001  format(' ERROR: atom',i6,' (',a4,a6,') has more than ',i2,
     -  ' neighbours and hydrogen bonds')
2002  format(' - increase the parameter MAXNEIG and recompile')
      end
