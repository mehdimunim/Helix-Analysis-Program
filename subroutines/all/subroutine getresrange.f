      subroutine getresrange(nsegm,indexs,isegno,ixres,iresno,ifres,
     -  label,llabel,numres,nslt,ires1save,ires2save,isegdef,iseg1,
     -  iseg2,listres,listlen,nrange,idefall,ionerange,maxrsd,maxrec,
     -  ihelp)
      dimension indexs(maxrec),isegno(maxrec),ixres(maxrec),
     -  iresno(maxrec),ifres(maxrsd),listres(maxrsd)
      character*(*) label
      common /logging/ logfile,ipredict
      character*80 line
c     print *,'GETRESRANGE nslt,numres,nsegm,ipredict=',
c    -  nslt,numres,nsegm,ipredict
      listlen=0
      nrange=0
      do while (.true.)
        if (nsegm .eq. 1 .and. ipredict .eq. 0) then
          call reverseindex(indexs,iresno,ifres,1,numres,maxrec)
9156      call getrange(ires1i,iresno(1),ires2i,
     -      iresno(nslt),incr,0,label,llabel,
     -      iresno(nslt),ihelp)
          ires1=indexs(ires1i)
          if (ires1 .eq. 0) write (6,2052) ires1i
          ires2=indexs(ires2i)
          if (ires2 .eq. 0) write (6,2052) ires2i
          if (ires1*ires2 .eq. 0) go to 9156
          iseg1=1
          iseg2=1
          print *,'ires1,ires2,iseg1,iseg2=',ires1,ires2,iseg1,iseg2
        else
          write (6,2039)
          line(1:28)='SEGMENT number of the first '
          lline=28
          line(lline+1:lline+llabel)=label(1:llabel)
          lline=lline+llabel
          if (isegdef .eq. 0) isegdef=1
9152      call getint(line,lline,isegdef,1,nsegm,iseg1,0)
          call findrange(isegno,1,nslt,iseg1,ifss,ilss,'segment',7,
     -      0,ifail)
          if (ifail .gt. 0) go to 9152
c         print *,'ifss,ilss=',ifss,ilss
c         print *,'iresno(ifss),iresno(ilss)=',iresno(ifss),iresno(ilss)
          write (6,2000) iseg1,iresno(ifss),iresno(ilss)
          line(1:7)='RESIDUE'
          call getint(line,lline,iresno(ifss),1,iresno(ilss),irn1,ihelp)
          line(1:28)='SEGMENT number of the last  '
          iseg2def=iseg1
          if (idefall .eq. 1) iseg2def=nsegm
9153      call getint(line,lline,iseg2def,1,numres,iseg2,ihelp)
          if (iseg2 .lt. iseg1) then
            print *,'ERROR: segment number of the last residue can not',
     -        ' be less than the first'
             go to 9153
          end if
          call findrange(isegno,1,nslt,iseg2,ifss,ilss,'segment',7,
     -      0,ifail)
          if (ifail .gt. 0) go to 9153
          write (6,2000) iseg2,iresno(ifss),iresno(ilss)
          line(1:7)='RESIDUE'
c         print *,'ifss,ilss=',ifss,ilss
c         print *,'iresno(ifss),iresno(ilss)=',iresno(ifss),iresno(ilss)
          call getint(line,lline,iresno(ilss),1,iresno(ilss),irn2,ihelp)
          call findsegres(isegno,iresno,ixres,1,nslt,iseg1,
     -      irn1,ia1,ires1,ifail1)
          call findsegres(isegno,iresno,ixres,ia1+1,nslt,
     -      iseg2,irn2,ia2,ires2,ifail2)
          if (ifail1+ifail2 .gt. 0) go to 9152
          if (ires2 .lt. ires1) then
            write (6,2001) ires2,ires1
            go to 9152
          end if
        end if
        do ir=ires1,ires2
          listlen=listlen+1
          listres(listlen)=ir
        end do
        nrange=nrange+1
        if (nrange .eq. 1) then
          ires1save=ires1
          ires2save=ires2
        else
          ires1save=min0(ires1,ires1save)
          ires2save=max0(ires2,ires2save)
        end if
        if (ionerange .eq. 1) return
        if (iseg1 .ne. iseg2) isegdef=0
        call askyn('Do you want to add an other range',33,1,-1,iadd,0,0)
        if (iadd .eq. 0) return
      end do
      return
2000  format(' Segment',i4,' includes residues ',i6,' to ',i6)
2001  format(' ERROR: residue index of the last residue (',i6,') is ',/,
     -  8x,'less than the first residue (',i6,')')
2039  format(' NOTE: segment numbers and residue ranges are printed ',
     -  'above')
2052  format(' ERROR: residue number',i6,' does not exist')
      end
