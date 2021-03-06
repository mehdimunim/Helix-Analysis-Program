      subroutine hbbridgeprint(nanchor,ianchor,lpath,nbridgetype,
     -  ibridgetype,maxbridgemem,line,index,iresno,inamcol1,inamcol2,
     -  irescol1,irescol2,iout,maxbhcount,maxhbtype,
     -  minbridgelenprint,minbridgepercprint,nframe,maxbridgelen,
     -  maxbridgetype,maxrec)
      dimension ianchor(nanchor),nbridgetype(maxbridgelen,nanchor),
     -  lpath(maxbridgetype,maxbridgelen,nanchor),iresno(maxrec),
     -  ibridgetype(maxbridgetype,maxbridgelen,nanchor),index(maxrec)
      character* 132 line(maxrec)
      dimension lensum(10)
c     print *,'HBBRIDGEPRINT maxbridgelen,maxbridgemem=',
c    -  maxbridgelen,maxbridgemem
      call zeroiti(lensum,0,maxbridgelen)
      if (minbridgelenprint .gt. 1) write (iout,1006) minbridgelenprint
      if (minbridgepercprint .gt. 0)
     -   write (iout,1007) minbridgepercprint,nframe
      write (iout,*)
      maxbrtyp=0
      maxbhcount=0
      maxhbtype=0
      do ia=1,nanchor
        ifrom=ianchor(ia)
        do id=1,maxbridgemem+1
          if (nbridgetype(id,ia) .gt. maxbrtyp)
     -       maxbrtyp=nbridgetype(id,ia)
          do in=1,nbridgetype(id,ia)
            ito=ibridgetype(in,id,ia)
            nbr=lpath(in,id,ia)
            if (id .ge. minbridgelenprint) then
              iperc=float(100*nbr)/float(nframe)
              if (iperc .ge. minbridgepercprint) then
                write (iout,1000) ifrom,iresno(ifrom),
     -            line(index(ifrom))(inamcol1:inamcol2),
     -            line(index(ifrom))(irescol1:irescol2),ito,iresno(ito),
     -            line(index(ito))(inamcol1:inamcol2),
     -            line(index(ito))(irescol1:irescol2),id,nbr
                if (id .lt. 1 .or. id .gt. 10) print *,
     -            'ia,ifrom,ito,nbr,id=',ia,ifrom,ito,nbr,id
                lensum(id)=lensum(id)+nbr
                if (id .eq. 1) then
                  maxhbtype=maxhbtype+1
                  if (nbr .gt. maxbhcount) maxbhcount=nbr
                end if
              end if
            end if
          end do
        end do
      end do
      write (iout,*)
      do id=1,maxbridgelen
        if (lensum(id) .gt. 0) write (iout,1001) id,lensum(id)
      end do
      write (6,1005) maxbrtyp,maxbridgetype
      write (iout,1005) maxbrtyp,maxbridgetype
      return
1000  format(' From (',i6,i6,')(',a,1x,a,') to (',i6,i6,
     -  ')(',a,1x,a,') of length ',i1,' #:',i5)
1001  format(' Number of paths of length',i2,'=',i9)
1005  format(/,' Maximum number of bridge types per anchor atom found=',
     -  i3,' limit=',i3)
1006  format(' Bridges with length < ',i1,' will not be printed')
1007  format(' Bridges that occur less than ',i3,'% of the',i6,
     -  ' frames will not be printed')
      end
