      subroutine savebond(i,j,iatnum,nneig,ineig,nhbneig,nhneig,
     -  ncneig,nnneig,npneig,nsneig,resnamj,atnamj,maxng,n,nosort,
     -  inamcol1,inamcol2,LEVTEST)
      dimension nneig(n),ineig(maxng,n),iatnum(n),nhbneig(n),
     -  nhneig(n),nnneig(n),ncneig(n),nsneig(n),npneig(n)
      character*4 atnamj
      character*8 resnamj
c     Save atom i as the neighbor of atom j
      if (LEVTEST .gt. 0) write (88,*) 'SAVEBOND i,j,=',i,j
      if (iatnum(j) .ge. 88 .and. iatnum(j) .le. 90 .or.
     -     iatnum(i) .lt. 88 .or. iatnum(i) .gt. 90) then
         maxngij=max0(nneig(i)+nhbneig(i),nneig(j)+nhbneig(j))
        if (maxngij .lt. maxng) then
          if (nneig(j) .gt. 0 .and. nosort .eq. 0) then
c           Keep the neighbor list sorted
            if (i .gt. ineig(nneig(j),j)) then
              ineig(nneig(j)+1,j)=i
            else
              nnj=nneig(j)
              do inn=1,nnj
                in=nnj-inn+1
                if (i .gt. ineig(in,j)) go to 11
              end do
              in=0
11            in0=in+1
              do inn=in0,nnj
                in=nnj-inn+in0
                ineig(in+1,j)=ineig(in,j)
              end do
              ineig(in0,j)=i
            end if
          else
            ineig(nneig(j)+1,j)=i
          end if
          nneig(j)=nneig(j)+1
          if (iatnum(i) .eq. 1) nhneig(j)=nhneig(j)+1
          if (iatnum(i) .eq. 6) ncneig(j)=ncneig(j)+1
          if (iatnum(i) .eq. 7) nnneig(j)=nnneig(j)+1
          if (iatnum(i) .eq. 16) nsneig(j)=nsneig(i)+1
          if (iatnum(i) .eq. 15) npneig(j)=npneig(i)+1
        else
          if (inamcol2 .ge. inamcol1) then
            write (6,2001) j,atnamj,resnamj,maxng
          else
            write (6,2000) j,maxng
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
