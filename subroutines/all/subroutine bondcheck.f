      subroutine bondcheck(iout,nfirst,n,iatnum,nneig,ineig,maxng,c,
     -  maxdist,line,irescol1,irescol2,inamcol1,inamcol2,index,rlim,
     -  nerr,maxrec)
      dimension nneig(n),ineig(maxng,n),c(3,n),iatnum(n),index(n),
     -  rlim(maxng)
      character* 132 line(maxrec)
      character*2 iatnm2
      common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
     -  mmatno(64),iatnm2(99)
c     Check for atoms exceeding their valence
      nerr=0
      nwarn=0
      nvalence=-1
      do i=nfirst,n
        if (iatnum(i) .gt. 0) nvalence=nval(iatnum(i))
        if (nneig(i) .gt. nvalence .and. nvalence .gt. 0) then
          do jj=1,nneig(i)
            j=ineig(jj,i)
            call decidebondcut(iatnum(i),iatnum(j),rlim(jj))
            rlim(jj)=sqrt(rlim(jj))
          end do
          write (iout,1000) i,iatnm2(iatnum(i)),
     -      line(index(i))(irescol1:irescol2),
     -      line(index(i))(inamcol1:inamcol2),nvalence,
     -      (ineig(j,i),line(index(ineig(j,i)))(irescol1:irescol2),
     -      line(index(ineig(j,i)))(inamcol1:inamcol2),
     -      sqrt(dist2(c(1,i),c(1,ineig(j,i)))),rlim(j),j=1,nneig(i))
          nerr=nerr+1
        end if
      end do
      if (nerr .gt. 0) then
        write (iout,1002) nerr
        write (6,1002) nerr
      else
        write (iout,*) 'No valence errors found'
      end if
c     Check for bonded atoms far apart
      do i=nfirst,n
        do jj=1,nneig(i)
          j=ineig(jj,i)
          if (iabs(i-j) .gt. maxdist .and. i .lt. j .and.
     -        (iatnum(i) .ne. 32 .or. iatnum(j) .ne. 32)) then
            write (iout,1001) i,line(index(i))(irescol1:irescol2),
     -        line(index(i))(inamcol1:inamcol2),
     -        j,line(index(j))(irescol1:irescol2),
     -        line(index(j))(inamcol1:inamcol2),maxdist
            nwarn=nwarn+1
          end if
        end do
      end do
      return
1000  format(' Atom',i7,1x,a2,' (',a4,1x,a4,') has more than',i2,
     -  ' bonds:',/,('   Bond to',i6,' (',a4,1x,a4,') d=',f4.2,
     -  ' A (threshold=',f4.2,' A)'))
1001  format(' Bond exists between atoms',i6,' (',a4,1x,a4,') and ',
     -  i6,' (',a4,1x,a4,')',/,' but they are more than ',i4,
     -  ' atoms apart')
1002  format(' !!! Number of valence errors found=',i3)
      end
