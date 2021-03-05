      subroutine findfg(n0,n,ian,nfgmem,ifgstr,ifgaix,indxfg,
     -  ixfg,ifgtyp,itypfg,nfg,iout,nneig,nneigh,nneiga,ineig,iwfg,
     -  inpcrdtyp,ioins,ic1,ic2,ir1,ir2,irn1,irn2,line,index,maxng,
     -  maxrec)
c#    MMC routine 164 lstmod: 07/11/96
c*****Assigns the appropriate fg types to the atoms in c
      dimension nneig(n),ineig(maxng,n),nneigh(n),nneiga(n),
     -  ian(n),nfgmem(n),ifgstr(n),ifgaix(n),ixfg(n),
     -  indxfg(n),ifgtyp(n),itypfg(n),index(n)
      character* 132 line(maxrec)
      character*2 iatnm2
      common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
     -  mmatno(64),iatnm2(99)
      character*1 sp
      character*4 namfcg
      common /connatdat/ ramax(99),ramax2(99),hlimfac,ianfg(99),
     -  namfcg(100),nrmw
      dimension in12(2)
      data sp /' '/,ityp /0/
c     c(1,i) : coordinates of atom i
c     n : number of atoms in c
c     icl(i) : type af atom i
c     ian(i) : atomic number of atom i
c     nfg: number of functional groups found
c     nfgmem(i): number of atoms in fcg i
c     ifgstr(i): start of member list for fcg i in indxat
c     itypfg(i): functional group type of atom i
c     ifgtyp(i): the functional group type if the i-th fcg before sorting
c     indxfg(i): the i-th atom belongs to the indxfg(i)-th fcg
c     ifgaix(i): the original index of the i-th atom in the atomlist
c     ixfg(ig) : the ig-th functional group after sorting
c     sorted by functional groups
c     nneig(i) : number of neighbours of atom i
c     ineig(1,i) : list of neighbours of atom i
c     nneigh(i) : number of hydrogen neighbours of atom i
      nfg=0
      nverr=0
      call zeroiti(itypfg,n0-1,n)
c     Check for "valence errors"
      nverr=0
      do i=n0,n
        if (nneiga(i) .gt. nval(ian(i)) .and. nval(ian(i)) .gt. 0) then
          if (inpcrdtyp .le. ioins) then
            write (iout,1001)
     -        i,line(index(i))(ic1:ic2),line(index(i))(ir1:ir2),
     -        line(index(i))(irn1:irn2),ian(i),nneiga(i),
     -        (ineig(j,i),line(index(ineig(j,i)))(ic1:ic2),
     -        line(index(ineig(j,i)))(ir1:ir2),
     -        line(index(ineig(j,i)))(irn1:irn2),j=1,nneiga(i))
          else
            write (iout,1002) i,ian(i),nneiga(i),
     -      (ineig(j,i),j=1,nneiga(i))
          end if
          itypfg(i)=100
          nfg=nfg+1
          indxfg(i)=nfg
          nverr=nverr+1
        end if
      end do
      if (nverr .gt. 0) then
        call askyn('Do you want to break any bond',29,1,1,ibreak,0,0)
        if (ibreak .gt. 0) then
          do while (.true.)
            call getintline(
     -        'Atomindices of the bond to be broken (0,0 to quit)',50,
     -        1,n,in12,2,00)
            if (in12(1) .eq. 0 .and. in12(2) .eq. 0) go to  100
            call breakbond(in12(1),in12(2),n0,n,nneig,ineig,nneiga,
     -        nneigh,ian,ifail,maxng)
            call breakbond(in12(2),in12(1),n0,n,nneig,ineig,nneiga,
     -        nneigh,ian,ifail,maxng)
          end do
        end if
      end if
c     if (nverr .gt. 0) go to 120
c     Search for O and P first
100   do i=n0,n
        if (ian(i) .eq. 8 .and. itypfg(i) .eq. 0 .and.
     -       nneiga(i) .gt. 0) then
c         Oxygen found
c          print *,'i,nneiga(i),ineig(1,i)=',i,nneiga(i),ineig(1,i)
          if (nneiga(i) .eq. 1 .and. ian(ineig(1,i)) .eq. 6) then
            inn1=ineig(1,i)
            nnn1=nneiga(inn1)
            noxy=0
            do j=1,nnn1
              if (ian(ineig(j,inn1)) .eq. 8 .and.
     -          nneiga(ineig(j,inn1)) .eq. 1) noxy=noxy+1
            end do
            if (noxy .eq. 2) then
c             COO- found typ=27
              nfg=nfg+1
              indxfg(inn1)=nfg
              itypfg(inn1)=27
              do j=1,nnn1
                in=ineig(j,inn1)
                if (ian(in) .eq. 8 .and. nneiga(in) .eq. 1 .and.
     -              itypfg(in) .eq. 0) then
                  indxfg(in)=nfg
                  itypfg(in)=27
                end if
              end do
            else
c             >C=O found,  typ=17
              nfg=nfg+1
              indxfg(i)=nfg
              indxfg(inn1)=nfg
              itypfg(i)=17
              itypfg(inn1)=17
            end if
          else if (nneiga(i) .eq. 2 .and. nneigh(i) .eq. 0) then
c           Ester oxygen found, typ=19 or 20 (for phospho ester)
            nfg=nfg+1
            indxfg(i)=nfg
            if (ian(ineig(1,i)) .ne. 15 .and. ian(ineig(2,i)) .ne. 15)
     -        itypfg(i)=19
            if (ian(ineig(1,i)) .eq. 15 .or. ian(ineig(2,i)) .eq. 15)
     -        itypfg(i)=20
          else if (nneiga(i) .eq. 2 .and. nneigh(i) .eq. 1) then
c           -OH found, typ=21
            nfg=nfg+1
            indxfg(i)=nfg
            itypfg(i)=21
            nn=nneiga(i)
            do j=1,nn
              in=ineig(j,i)
              if (ian(in) .eq. 1 .and. itypfg(in) .eq. 0) then
                indxfg(ineig(j,i))=nfg
                itypfg(ineig(j,i))=21
              end if
            end do
          end if
        else if (ian(i) .eq. 15 .and. itypfg(i) .eq. 0) then
c         Phosphorus found  (>PO2 : ityp 22)
          nfg=nfg+1
          indxfg(i)=nfg
          itypfg(i)=22
          iong=0
          nn=nneiga(i)
          do j=1,nn
            in=ineig(j,i)
            if (ian(in) .eq. 8 .and. itypfg(in) .eq. 0 .and.
     -        nneiga(in) .eq. 1) then
              iong=iong+1
              indxfg(in)=nfg
              itypfg(in)=22
            end if
          end do
          if (iong .ne. 2) write (iout,1000) i,iong
        end if
      end do
c     Search for carbon and nitrogen next
      do i=n0,n
        if (ian(i) .eq. 6 .and. itypfg(i) .eq. 0) then
c         Unassigned carbon found
          nfg=nfg+1
          ityp=(10-(nneiga(i)*(nneiga(i)+1))/2)+nneigh(i)+1
          if (ityp .lt. 0) ityp=99
          indxfg(i)=nfg
          itypfg(i)=ityp
          nn=nneiga(i)
          if (nn .gt. 0) then
            do j=1,nn
              in=ineig(j,i)
              if (ian(in) .eq. 1 .and. itypfg(in) .eq. 0) then
                indxfg(in)=nfg
                itypfg(in)=ityp
              end if
            end do
          end if
        else if (ian(i) .eq. 7 .and. itypfg(i) .eq. 0) then
c         Unassigned nitrogen found
          nfg=nfg+1
          if (nneiga(i) .lt. 4)
     -      ityp=(6-(nneiga(i)*(nneiga(i)+1))/2)+nneigh(i)+11
c         Four neighbours, assumed to be positively charged
          if (nneiga(i) .eq. 4) ityp=23+nneigh(i)
          indxfg(i)=nfg
          itypfg(i)=ityp
          nn=nneiga(i)
          if (nn .gt. 0) then
            do j=1,nn
              in=ineig(j,i)
              if (ian(in) .eq. 1 .and. itypfg(in) .eq. 0) then
                indxfg(in)=nfg
                itypfg(in)=ityp
              end if
            end do
          end if
        end if
      end do
c     Label hydrogens on carbonyl
      do i=n0,n
        if (itypfg(i) .eq. 0 .and. nneiga(i) .gt. 0) then
c         Hydrogen on a C=O is type 18
          if (ian(i) .eq. 1 .and. itypfg(ineig(1,i)) .eq. 17) then
            nfg=nfg+1
            indxfg(i)=nfg
            itypfg(i)=18
          end if
        end if
      end do
c     Search for -S- (type 28) and -SH (type 29)
      do i=n0,n
        if (itypfg(i) .eq. 0 .and. ian(i) .eq. 16) then
c         Sulphur found
          if (nneiga(i) .eq. 2) then
c           -S-, -SH or HSH found
            if (nneigh(i) .gt. 0) then
c             -SH or HSH found
              itp=29
              nfg=nfg+1
              indxfg(i)=nfg
              itypfg(i)=itp
              do ing=1,2
                in=ineig(ing,i)
                if (ian(in) .eq. 1 .and. itypfg(in) .eq. 0) then
                  indxfg(in)=nfg
                  itypfg(in)=itp
                end if
              end do
            else
c             -S- found
              nfg=nfg+1
              indxfg(i)=nfg
              itypfg(i)=28
            end if
          end if
        end if
      end do
c     Assign the single atom funcional groups
      do i=n0,n
        if (itypfg(i) .eq. 0) then
          ianfgi=ianfg(ian(i))
          if (ianfgi .gt. 0) then
c           If no neighbours, must be an ion
            if (nneiga(i) .eq. 0) ianfgi=ianfgi+1
            nfg=nfg+1
            indxfg(i)=nfg
            itypfg(i)=ianfgi
          end if
        end if
      end do
c     Label all unassigned atoms
      do i=n0,n
        if (itypfg(i) .eq. 0) then
          nfg=nfg+1
          indxfg(i)=nfg
          itypfg(i)=99
        end if
      end do
c     Sort atoms by functional groups
      nmem=0
      do if=1,nfg
        ifgstr(if)=nmem+1
        do ia=n0,n
          if (indxfg(ia) .eq. if) then
            nmem=nmem+1
            ifgaix(nmem)=ia
          end if
          nfgmem(if)=nmem-ifgstr(if)+1
        end do
        ifgtyp(if)=itypfg(ifgaix(ifgstr(if)))
      end do
c     index (sort) fcg's by type
      nfgn=0
      do it=1,100
        do if=1,nfg
          if (ifgtyp(if) .eq. it) then
            nfgn=nfgn+1
            ixfg(nfgn)=if
          end if
        end do
      end do
c      if (nfg .gt. 1) then
cc       Sort functional groups
c        do 50 i=1,nfg
c          j1=i+1
c          do 50 j=j1,nfg
c            if (ifgtyp(i) .gt. ifgtyp(j)) then
c              ii=ifgtyp(i)
c              ifgtyp(i)=ifgtyp(j)
c              ifgtyp(j)=ii
c              ii=ifgstr(i)
c              ifgstr(i)=ifgstr(j)
c              ifgstr(j)=ii
c              ii=nfgmem(i)
c              nfgmem(i)=nfgmem(j)
c              nfgmem(j)=ii
c            end if
c50      continue
c      end if
      if (iwfg .gt. 0) then
c       Print list
        do igg=1,nfg
          ig=ixfg(igg)
          ig1=ifgstr(ig)
          ig2=ig1+nfgmem(ig)-1
          do iat=ig1,ig2
            iaa=ifgaix(iat)
            nna=nneig(iaa)
            if (inpcrdtyp .le. ioins) then
              write (iwfg,2045) igg,namfcg(ifgtyp(ig)),
     -          iaa,line(index(iaa))(ic1:ic2),
     -          line(index(iaa))(ir1:ir2),
     -          line(index(iaa))(irn1:irn2),(ineig(in,iaa),
     -          line(index(ineig(in,iaa)))(ic1:ic2),
     -          line(index(ineig(in,iaa)))(ir1:ir2),
     -          line(index(ineig(in,iaa)))(irn1:irn2),in=1,nna)
            else
              write (iwfg,2044) igg,namfcg(ifgtyp(ig)),
     -          iaa,iatnm2(ian(iaa)),(sp,ineig(in,iaa),
     -          iatnm2(ian(ineig(in,iaa))),in=1,nna)
            end if
          end do
        end do
      end if
      return
1000  format(' ERROR: phosporus(',i4,') has',i2,' oxygens')
1001  format(' ERROR: atom',i6,' (',a,1x,a,1x,a,' atomic no=',i2,')',
     -    ' has',i3,' neighbours: ',/,3(i6,' (',a,1x,a,1x,a,') '))
1002  format(' ERROR: atom',i6,', atomic no=',i2,
     -    ' has',i3,' neighbours: ',(10i6))
2044  format(' Fg',i6,2x,a4,' (',i5,1x,a2,')',' neighbours:',a1,
     -  '(',i5,1x,a2,')',a1,'(',i5,1x,a2,')',a1,'(',i5,1x,a2,')',a1,
     -  /,(5x,5('(',i5,1x,a2,')',a1)))
2045  format(' Fg',i6,2x,a4,' (',i5,1x,a4,1x,a4,1x,a,') Neigbours:',/,
     -  (9x,3('(',i5,1x,a4,1x,a4,1x,a,') ')))
      end
