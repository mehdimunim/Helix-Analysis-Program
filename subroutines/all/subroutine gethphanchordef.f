      subroutine gethphanchordef(line,index,nslt,iresno,iatnum,charge,
     -  indexa,indexn,indexo,nneig,ineig,nhneig,nhneigmin,nanchor,
     -  ianchor,ianchor2,iselfanc,nosameseg,iallheavy,iiq1,iiq2,
     -  inpcrdtyp,iobpdb,iocpdb,icharges,ibondtype,iw0,segid4,
     -  molsltlim,nsegslt,bondname,lbondname,ifail,maxneig,maxrec,
     -  maxanchorlist)
      dimension index(maxrec),charge(maxrec),iatnum(maxrec),
     -  indexa(maxrec),indexn(maxrec),indexo(maxrec),iresno(maxrec),
     -  ianchor(maxanchorlist),nneig(maxrec),ineig(maxneig,maxrec),
     -  nhneig(maxrec),molsltlim(3,nsegslt)
      character*4 segid4(nsegslt)
      character* 132 line(maxrec)
      character*1 ansrun
      character*(*) bondname
c     Establish solute atoms that are hydrophobic bond  anchors
c     print *,'GETHPHANCHOR maxrec,maxanchorlist=',
c    -         maxrec,maxanchorlist
      nanchor=0
      iqsel=0
      qmin=0.0
      call zeroiti(indexa,0,nslt)
      nanchorr=0
9134  ansrun=' '
      iall=0
      do while (iall .eq. 0 .and. ansrun .ne. 'q')
        call quiz(ansrun,iansrun,' ',' ',1,
     -    'hydrophobic/salt bridge/contact anchor atoms',44,0,5,6,75)
        ifail=0
        call definelist(ansrun,nslt,nanchorr,indexn,indexo,nsegslt,
     -    segid4,iresno,molsltlim,bondname,lbondname,iall,maxanchorlist)
      end do
      iqfsel=0
      if (icharges .gt. 0)
     -  call askyn('Do you want to filter the list by charge',40,
     -    1,-1,iqfsel,0,0)
      if (iqfsel .gt. 0) then
        if (icharges .eq. 1) then
          if (ischarmm(inpcrdtyp) .eq. 1)
     -      print *,'Charges are read from the WEIGHT column'
          if (inpcrdtyp .eq. iobpdb .or. inpcrdtyp .eq. iocpdb) print *,
     -      'Charges are read from the TEMPERATURE FACTOR column'
        end if
        call getreal('Maximum (absolute) charge for an anchor atom',
     -    44,999999.0,qmax,1,0)
        iall=0
      end if
c     indexa will be nonzero for possible anchor atoms, negative if not anchor
      if (ibondtype .eq. 2) then
c       Hydrophobic 'bond'
8010    nhpc=0
        do ia=1,nslt
          if (iatnum(ia) .eq. 6 .and. nhneig(ia) .ge. nhneigmin) then
            indexa(ia)=-1
            nhpc=nhpc+1
          end if
        end do
        if (nhneigmin .gt. 0 .and. nhpc .eq. 0) then
          print *,'No H bonded to C was found'
          call askyn('Do you want to just use all carbons',35,1,+1,i0,0,
     -      0)
          if (i0 .eq. 1) then
            nhneigmin=0
            go to 8010
          else
            ifail=1
            return
          end if
        end if
        icofil=0
        if (iallheavy .eq. 0)
     -    call askyn('Do you want to filter out >C=O and C-OH carbons',
     -      47,1,+1,icofil,0,0)
        if (icofil .gt. 0) then
          do ia=1,nslt
            if (indexa(ia) .lt. 0) then
              idrop=0
              do in=1,nneig(ia)
                inc=ineig(in,ia)
                if (iatnum(inc) .eq. 8) then
                  if (nneig(inc) .eq. 1 .or. nhneig(inc) .gt. 0) idrop=1
                end if
              end do
              if (idrop .eq. 1) indexa(ia)=0
            end if
          end do
        end if
      else
c       Heavy-atom contact
        do ia=1,nslt
          if (iatnum(ia) .ne. 1) indexa(ia)=-1
        end do
      end if
      do ja=1,nanchorr
        ia=indexn(ja)
c       print *,'Anchor #',ja,'=',ia,' indexa(ia)=',indexa(ia),
c    -    ' iatnum(ia)=',iatnum(ia)
        if (indexa(ia) .ne. 0) then
          iasel=1
          if (iqfsel .gt. 0) then
            if (icharges .eq. 1) then
              call readreal(line(index(ia)),iiq1,iiq2,qa)
            else
              qa=charge(ia)
            end if
            if (qa .lt. -qmax .or. qa .gt. qmax) iasel=0
          end if
          if (iasel .gt. 0) then
            indexa(ia)=-indexa(ia)
            if (nanchor .eq. maxanchorlist) then
              print *,'Anchor list limit (',maxanchorlist,
     -          ') is reached'
              nanchorr=ja
            else
              nanchor=nanchor+1
              ianchor(nanchor)=ia
            end if
          end if
        else if (ansrun .eq. 'l' .and. iallheavy .eq. 0) then
          print *,'NOTE: atom ',ia,' is not a hydrophobic carbon'
        end if
      end do
      if (nanchor .eq. 0) then
        print *,'PROBLEM: No anchor atoms were found'
        call askyn('Do you want to repeat the selection',35,
     -    1,+1,iselrep,0,0)
        if (iselrep .eq. 1) go to 9134
        ifail=1
        return
      end if
      call getanchormod(ianchor2,iselfanc,nosameseg,iall,iw0)
      do ia=1,nanchor
        indexa(ianchor(ia))=1
      end do
      return
      end
