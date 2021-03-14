      subroutine hbbridge(nanchor,ianchor,indexa,ianchor2,iselfanc,
     -  lpath,nbridgetype,ibridgetype,maxbridgemem,n,nhbneig,ineig,list,
     -  nnblist,iparent,ixres,resnames,brslv,nabr,nrescol,nmc,ifail,
     -  listbridge,iout,maxng,maxbridgelen,maxbridgetype,
     -  minbridgelenprint,maxrsd,maxrec)
      character*8 resnames
      dimension ianchor(nanchor),indexa(n),nhbneig(n),ineig(maxng,n),
     -  lpath(maxbridgetype,maxbridgelen,nanchor),
     -  nbridgetype(maxbridgelen,nanchor),
     -  ibridgetype(maxbridgetype,maxbridgelen,nanchor),resnames(maxrsd)
      dimension list(maxrec),nnblist(maxrec),iparent(maxrec),
     -  ixres(maxrec),ip(10),iflist(10),illist(10)
      character*8 brslv
      character*80 line
c     write (iout,*) 'ianchor2,iselfanc,n,ixres(1),ixres(n)=',
c    -  ianchor2,iselfanc,n,ixres(1),ixres(n)
c     print *,'HBBRIDGE nanchor=',nanchor,' maxrec=',maxrec
      ifail=0
      do iia=1,nanchor
        do ia=1,ixres(n)
          iparent(ia)=-1
        end do
        iat=ianchor(iia)
c       write (iout,2311) iat,indexa(iat),resnames(ixres(iat)),
c    -    (ineig(maxng+1-ia,iat),ia=1,nhbneig(iat))
c2311    format(' HBBR ianch=',i6,' indexa=',i4,' resn=',a,' ihb=',5i6)
        if (indexa(iat) .eq. 0) then
          print *,'PROGRAM ERROR: ',iia,'th anchorlist=',iat,
     -      ' is not marked as anchor'
          stop
        end if
        list(1)=iat
        ll=1
        ixp=ixres(iat)
c       Grow tree at depth maxbrlen (4)
        ll1=1
        iflist(1)=ll1
        illist(1)=ll
        do it=1,maxbridgemem+1
          nn=0
          do il=ll1,ll
            ia=list(il)
c           ixp is the parent residue (solvent molecule) number
            if (ia .ne. iat) ixp=ixres(ia)
            if (ixp .le. 0) then
              write (iout,1005) iat,ia,it,il,ixp
              write (6,1005) iat,ia,it,il,ixp
              stop
            end if
            do in=1,nhbneig(ia)
              ian=ineig(maxng-in+1,ia)
c             Drop if it is a loop back or already in the current list
              do j=1,illist(min0(2,it))
                if (list(j) .eq. ian) go to 100
              end do
c             list is only sorted at each level
              do itt=3,it
c               print *,'it,itt,iflist(itt),illist(itt)=',
c    -             it,itt,iflist(itt),illist(itt)
                call findixsort(list,iflist(itt),illist(itt),ian,
     -            ixian,itryl)
                if (ixian .gt. 0) go to 100
              end do
c             do j=1,il-1
c               if (nnblist(j) .eq. ian) go to 100
c             end do
              call findixsort(nnblist,1,nn,ian,ixian,itrynnb)
              if (ixian .gt. 0) go to 100
c             If anchor atom then finish bridge and gather statistics
              if (resnames(ixres(ian))(1:nrescol) .ne.
     -            brslv(1:nrescol)) then
                idoit=1
c               See if end is also an anchor
                if (ianchor2 .eq. 1) then
                  if (indexa(ian) .lt. 1) idoit=0
                else
                  if (iselfanc .eq. 0 .and. indexa(ian) .eq. 1) idoit=0
                end if
                if (idoit .eq. 1 .and.
     -              (indexa(ian) .lt. 1 .or. ian .gt. iat)) then
                  if (nbridgetype(it,iia) .eq. 0) then
                    nbridgetype(it,iia)=1
                    ibridgetype(1,it,iia)=ian
                    iend=1
                    go to 300
                  end if
                  call findixsort(ibridgetype(1,it,iia),1,
     -               nbridgetype(it,iia),ian,ixian,itry)
                  if (ixian .gt. 0) then
                    iend=ixian
                    go to 300
                  end if
c                 do ial=1,nbridgetype(it,iia)
c                   if (ibridgetype(ial,it,iia) .eq. ian) then
c                     iend=ial
c                     go to 300
c                   end if
c                 end do
                  if (nbridgetype(it,iia) .lt. maxbridgetype) then
c                   Not found, add to list
c                   Shift ibridgetype,lpath to maintain them sorted
                    ixian=itry
                    do ib=nbridgetype(it,iia),ixian,-1
                      ibridgetype(ib+1,it,iia)=ibridgetype(ib,it,iia)
                      lpath(ib+1,it,iia)=lpath(ib,it,iia)
                    end do
                    nbridgetype(it,iia)=nbridgetype(it,iia)+1
                    ibridgetype(ixian,it,iia)=ian
                    iend=ixian
c                   ibridgetype(nbridgetype(it,iia),it,iia)=ian
c                   iend=nbridgetype(it,iia)
                  else
                    write (6,1006) maxbridgetype
                    write (iout,1006) maxbridgetype
                    if (nmc .gt. 0) then
                      write (6,1008) nmc
                      write (iout,1008) nmc
                      ifail=1
                      return
                    else
                      percdone=100.0*float(iia)/float(nanchor)
                      write (6,1009) percdone
                      write (iout,1009) percdone
                      stop
                    end if
                  end if
300               lpath(iend,it,iia)=lpath(iend,it,iia)+1
                  if (listbridge .gt. 0) then
c                   write (iout,1002) ibridgetype(iend,it,iia)
                    ixptb=ixp
                    ip(it)=ixptb
                    iit=it
                    do while (iit .gt. 2)
                      iit=iit-1
c                     write (iout,*) 'TBCK it,iit,ixp,iparent(ixp)=',
c    -                 it,iit,ixptb,iparent(ixptb)
                      ixptb=iparent(ixptb)
                      ip(iit)=ixptb
                      if (ixp .eq. -1) then
                        write (6,1007) iat,it
                        write (iout,1007) iat,it
                        stop
                      end if
                    end do
                    if (nmc .eq. 0) then
                      ic0=0
                    else
                      ic0=10
                      write (line(1:10),1004) nmc
                    end if
                    write (line(ic0+1:ic0+18),1000) iat
                    ic=ic0+18
                    if (it .gt. minbridgelenprint) then
                      do iit=2,it
                        write (line(ic+1:ic+17),1001) ip(iit)
                        ic=ic+17
                        if (iit .eq. 2 .and. iit .lt. it) then
                          write (iout,1003) line(1:ic)
                          ic=ic0+18
                          call blankout(line,1,ic)
                        end if
                      end do
                    end if
                    write (line(ic+1:ic+15),1002)
     -                ibridgetype(iend,it,iia)
                    write (iout,1003) line(1:ic+15)
                  end if
c                  write (77,8877) il,ia,iend,it,iia,lpath(iend,it,iia)
c8877              format(' il,ia,iend,it,iia=',5i4,' lpath=',i3)
                end if
              else
c               New atom, add to list of next level all atoms of this solvent
                isolv=ixres(ian)
                iparent(isolv)=ixp
c               write (iout,*) 'ADD ian,isolv,iparent(isolv)=',
c    -                              ian,isolv,iparent(isolv)
c               Add all atoms of this residue to nnblist
                iaa=ian
                do while (iaa .gt. 1 .and. ixres(iaa) .eq. ixres(ian))
                  iaa=iaa-1
                end do
                if (iaa .eq. 1 .and. ixres(iaa) .eq. ixres(ian))
     -            iaa=iaa-1
                ifr=iaa+1
                iaa=ian
                do while (iaa .lt. n .and. ixres(iaa) .eq. ixres(ian))
                  iaa=iaa+1
                end do
                if (iaa .eq. n .and. ixres(iaa) .eq. ixres(ian))
     -            iaa=iaa+1
                ilr=iaa-1
                if (ilr-ifr+1 .ne. nabr) then
                  write (6,1010) brslv,ixres(ian),ilr-ifr+1,nabr
                  write (iout,1010) brslv,ixres(ian),ilr-ifr+1,nabr
                  stop
                end if
c               write (iout,*) 'Adding itrynnb,ifr,ilr,nn=',
c    -                          itrynnb,ifr,ilr,nn
                if (nn .gt. 0) then
                  do is=nn,itrynnb,-1
                    nnblist(is+nabr)=nnblist(is)
                  end do
                end if
                do isv=1,nabr
                  nnblist(itrynnb-1+isv)=ifr-1+isv
                end do
                nn=nn+nabr
c                 do isv=ifr,ilr
c                   nnblist(nn+isv-ifr+1)=isv
c                 end do
c               nn=nn+ilr-ifr+1
              end if
100           continue
            end do
          end do
c          if (nn .gt. 0) write (iout,6633) (nnblist(kk),kk=1,nn)
c6633      format(' nnblist=',10i6)
          call trnsfi(list(ll+1),nnblist,nn)
          ll1=ll+1
          ll=ll+nn
          iflist(it+1)=ll1
          illist(it+1)=ll
        end do
      end do
      return
1000  format(' Anchor 1:',i5,' - ')
1001  format(' solvent',i6,' - ')
1002  format(' Anchor 2:',i5)
1003  format(a)
1004  format(i10)
1005  format(' PROGRAM ERROR at iat,ia,it,il=',2i5,i3,i2,
     -  ': ixres(ia)=',i4)
1006  format(' ERROR: maximum number of bridge end types (',i4,
     -  ') is exceeded',/,8x,'Reduce the number of bridge members or',/,
     -  8x,'increase MAXBRIDGETYPE')
1007  format(' PROGRAM ERROR in traceback iat=',i4,' it=',i1)
1008  format(8x,'Last frame processed=',i9)
1009  format(8x,'Processed ',f5.1,' % of the anchors')
1010  format(' ERROR: bridge residue ',a,' #',i5,' has',i4,' atoms ',
     -  'instead of',i3,/,' - check the input STRUCTURE file')
      end
