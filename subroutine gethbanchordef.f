      subroutine gethbanchordef(line,index,nslt,ixres,iresno,iatnum,
     -  indexa,indexo,nanchor,ianchor,iiq1,iiq2,inpcrdtyp,iobpdb,iocpdb,
     -  icharges,qmin,iw0,inamcol1,resnames,brslv,nrescol,segid4,
     -  molsltlim,nsegslt,isc,nneig,ineig,label,llabel,ianchor2,
     -  iselfanc,nosameseg,iqfsel2,iquizhelp,ifail,maxng,maxrsd,maxrec,
     -  maxanchorlist)
      dimension index(maxrec),indexa(maxrec),
     -  indexo(maxrec),ixres(maxrec),iresno(maxrec),iatnum(maxrec),
     -  isc(maxrec),ianchor(maxanchorlist),molsltlim(3,maxrsd),
     -  nneig(maxrec),ineig(maxng,maxrec)
      character*4 segid4(nsegslt)
      character*8 resnames(maxrsd)
      character* 132 line(maxrec)
      character*(*) brslv,label
      character*4 atnam,bridgeats(200)
      character*1 ansrun
      character*4 bbhbats(6)
      data bbhbats /'O   ',' O  ','H   ',' H  ','HN  ',' HN '/,
     -  nbbhbats /6/
c     print *,'GETHBANCHORDEF maxanchorlist=',maxanchorlist
c     Establish solute atoms that are H-bond (chain) anchors
      nanchor=0
      iqsel=0
      iallsel=0
      qmin=0.0
      iskipbb=0
      iselbb=0
      nanchorr=0
9134  ansrun=' '
      iall=0
      do while (iall .eq. 0 .and. ansrun .ne. 'q')
        call quiz(ansrun,iansrun,' ',label,llabel,
     -    'hydrogen bond anchor atoms',26,0,5,6,iquizhelp)
        ifail=0
        if (ansrun .eq. 'b' .or. ansrun .eq. 't') then
          if (ansrun .eq. 'b') then
c           Select backbone atoms
            do i=1,nbbhbats
              bridgeats(i)=bbhbats(i)
            end do
            nbridgeats=nbbhbats
            iselbb=1
          else
c           Input atom name list
            call getnamelist(bridgeats,4,nbridgeats,'Anchor atom',11,
     -        200)
          end if
          write (6,7731) label(1:llabel),(bridgeats(i),i=1,nbridgeats)
          do ia=1,nslt
            if (resnames(ixres(ia))(1:nrescol) .ne.
     -          brslv(1:nrescol)) then
              atnam=line(index(ia))(inamcol1:inamcol1+3)
c             print *,'ia,atnam=',ia,atnam
              do it=1,nbridgeats
                if (atnam .eq. bridgeats(it)) then
c                 Anchor atom found
                  if (nanchor .lt. maxanchorlist) then
                    nanchor=nanchor+1
                    ianchor(nanchor)=ia
                    go to 9133
                  else
                    write (6,2144) maxanchorlist,label(1:llabel)
                    go to 9132
                  end if
                end if
              end do
9133          continue
            end if
          end do
        else if (ansrun .eq. 'c') then
c         Select by partial charge
          iqsel=1
        else
          call definelist(ansrun,nslt,nanchor,ianchor,indexo,nsegslt,
     -      segid4,iresno,molsltlim,label,llabel,iallsel,maxanchorlist)
        end if
      end do
      if (iselbb .eq. 0)
     -  call askyn('Do you want to omit protein backbone atoms',42,
     -    1,-1,iskipbb,0,0)
      iqfsel=0
      if (iqsel .eq. 0 .and. icharges .gt. 0) then
        call askyn('Do you want to filter the list by charge',40,
     -    1,-1,iqfsel,0,0)
      end if
      if (iskipbb+iqfsel .gt. 0) iallsel=0
      if (iqsel+iqfsel .gt. 0) then
        if (iiq1 .gt. iiq2) then
          print *,'This input format has no charge information'
          go to 9134
        end if
        if (ischarmm(inpcrdtyp) .eq. 1)
     -    print *,'Charges are read from the WEIGHT column'
        if (inpcrdtyp .eq. iobpdb .or. inpcrdtyp .eq. iocpdb) print *,
     -    'Charges are read from the TEMPERATURE FACTOR column'
        call getreal('Minimum (absolute) charge for an anchor atom',44,
     -    999999.0,qmin,1,0)
        call zeroiti(indexa,0,nslt)
        nzrq=0
        qsum=0.0
        do ia=1,nslt
          if (resnames(ixres(ia))(1:nrescol) .ne.
     -        brslv(1:nrescol)) then
            call readreal(line(index(ia)),iiq1,iiq2,qa)
            if (qa .eq. 0.0) nzrq=nzrq+1
            qsum=qsum+qa
            if (qa .lt. -qmin .or.
     -        (qa .gt. qmin .and. iatnum(ia) .eq. 1)) indexa(ia)=-1
          end if
        end do
        if (nzrq .eq. nslt) then
          print *,'Input structure had all zero charges'
          write (6,2145)
          go to 9134
        else
          print *,'Sum of solute charges read=',qsum
          if (nslt .gt. 10 .and. abs(qsum) .gt. alog(float(nslt)))
     -      write (6,2145)
        end if
        if (iqsel .eq. 1) then
c         Select from all solute atoms
          do ia=1,nslt
            if (indexa(ia) .eq. -1) then
              if (nanchor .lt. maxanchorlist) then
                nanchor=nanchor+1
                ianchor(nanchor)=ia
              else
                write (6,2144) maxanchorlist,label(1:llabel),' ',ia
                go to 9132
              end if
            end if
          end do
        else
c         Reduce list
          do iaa=1,nanchor
            if (indexa(ianchor(iaa)) .eq. 0) ianchor(iaa)=0
          end do
          ndel=0
          do ia=1,nanchor
            if (ianchor(ia) .eq. 0) then
              ndel=ndel+1
            else
              ianchor(ia-ndel)=ianchor(ia)
            end if
          end do
          nanchor=nanchor-ndel
          print *,'Number of anchor atoms filtered out=',ndel
          if (nanchor .eq. 0 .and. ndel .gt. 0) then
            print *,'PROBLEM: No anchor atoms were left after filtering'
            call askyn('Do you want to repeat the selection',35,
     -        1,+1,iselrep,0,0)
            if (iselrep .eq. 1) go to 9134
            ifail=1
            return
          end if
        end if
      end if
      call getanchormod(ianchor2,iselfanc,nosameseg,iallsel,iw0)
      if (qmin .ne. 0.0) write (iw0,*) 'Anchor atoms were filtered ',
     -  'by a charge threshold of ',qmin
      iqfsel2=0
      if (iqsel+iqfsel .gt. 0)
     -  call askyn(
     -    'Do you want to apply the charge filter to all atoms',51,
     -    1,+1,iqfsel2,0,0)
      if (iqfsel2 .eq. 1) write (iw0,*)
     -  'Charge filter will be applied to all putative anchor atoms'
      if (nanchor .gt. maxanchorlist) then
        write (6,2144) maxanchorlist,label(1:llabel),' ',ia
      end if
      ndelbb=0
      if (iskipbb .gt. 0) then
        do ia=1,nanchor
          ib=ianchor(ia)
          if (isc(ib) .eq. 0) then
c           BB atom - skip
            ndelbb=ndelbb+1
          else
             ianchor(ia-ndelbb)=ianchor(ia)
          end if
        end do
        if (ndelbb .gt. 0) then
          nanchor=nanchor-ndelbb
          write (6,2147) 'backbone atoms',ndelbb
          write (iw0,2147) 'backbone atoms',ndelbb
        end if
      end if
      write (iw0 ,*) 'ANCHOR LIST before ALIPH removal'
      call condenselist(ianchor,nanchor,0,iw0 )
c     Remove carbons and aliphatic hydrogens
      ndelhc=0
      do ia=1,nanchor
        if (nneig(ia) .eq. 0) then
          ianchor(ia-ndelhc)=ianchor(ia)
        else
          iaa=ianchor(ia)
          if (iatnum(iaa) .eq. 6 .or.
     -      (iatnum(iaa) .eq. 1 .and. iatnum(ineig(1,iaa)) .eq. 6)) then
            ndelhc=ndelhc+1
          else
            ianchor(ia-ndelhc)=ianchor(ia)
          end if
        end if
      end do
      if (ndelhc .gt. 0) then
        nanchor=nanchor-ndelhc
        write (6,2147) 'backbone atoms',ndelhc
        write (iw0,2147) 'carbons and aliphatic hydrogens',ndelhc
      end if
      if (ndelbb+ndelhc .gt. 0) write (iw0,2148) nanchor
      do ia=1,nanchor
        indexa(ianchor(ia))=1
      end do
      return
9132  if (nanchor .eq. 0) then
        print *,'PROBLEM: No anchor atoms were found'
        call askyn('Do you want to repeat the selection',35,
     -    1,+1,iselrep,0,0)
        if (iselrep .eq. 1) go to 9134
        ifail=1
        return
      end if
      return
2144  format(' Progam is limited to ',i5,1x,a,' anchor ',
     -  'atoms.',/,' Redimension or break up the run ',a,/,
     -  ' Last solute atom used is the ',i6,'-th')
2145  format(' Most likely, the charge field does not contain the ',
     -  'charges')
2147  format(' Number of ',a,' filtered out=',i5)
2148  format(' Number of hydrogen-bonding anchor atom left=',i5)
7731  format(1x,a,' anchor atom names:',/,(10(a4,1x)))
      end
