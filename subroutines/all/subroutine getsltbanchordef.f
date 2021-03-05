      subroutine getsltbanchordef(line,index,nslt,iresno,iatnum,charge,
     -  indexa,indexn,indexo,chargesum,nhneig,ineig,nanchor,ianchor,
     -  ianchor2,iselfanc,nosameseg,iiq1,iiq2,inpcrdtyp,iobpdb,iocpdb,
     -  icharges,iw0,irescol1,irescol2,inamcol1,inamcol2,segid4,
     -  molsltlim,nsegslt,isegno,ifail,maxneig,maxrec,maxanchorlist)
      dimension index(maxrec),charge(maxrec),indexa(maxrec),
     -  indexn(maxrec),indexo(maxrec),iresno(maxrec),chargesum(maxrec),
     -  iatnum(maxrec),ianchor(maxanchorlist),nhneig(maxrec),
     -  ineig(maxneig,maxrec),molsltlim(3,nsegslt),isegno(maxrec)
      character*4 segid4(nsegslt)
      character* 132 line(maxrec)
      character*1 ansrun
      character*8 saltatname(100),atnam
c     Establish solute atoms that are hydrophobic bond  anchors
c     print *,'GETHPHANCHOR maxneig,maxrec,maxanchorlist=',
c    -         maxneig,maxrec,maxanchorlist
      nsaltats=0
      qmin=0.0
      call zeroiti(indexa,0,nslt)
      call zeroit(chargesum,nslt)
      call quiz(ansrun,iansrun,' ',' ',1,
     -  'salt-bridge definition mode',27,0,5,6,0)
      if (ansrun .eq. 'c' .or. ansrun .eq. 'b') then
c       Atom selection by charge
        if (iiq1 .gt. iiq2) then
          print *,'This input format has no charge information'
        else
          if (icharges .eq. 1) then
            if (ischarmm(inpcrdtyp) .eq. 1)
     -        print *,'Charges are read from the WEIGHT column'
            if (inpcrdtyp .eq. iobpdb .or. inpcrdtyp .eq. iocpdb)
     -        print *,
     -          'Charges are read from the TEMPERATURE FACTOR column'
          end if
          call getreal('Minimum (absolute) charge for an anchor atom',
     -      44,999999.0,qmin,1,0)
          do ia=1,nslt
            if (icharges .eq. 1) then
              call readreal(line(index(ia)),iiq1,iiq2,chargesum(ia))
            else
              chargesum(ia)=charge(ia)
            end if
          end do
          do ia=1,nslt
            if (iatnum(ia) .ne. 1 .and. iatnum(ia) .ne. 6) then
              do in=1,nhneig(ia)
                chargesum(ia)=chargesum(ia)+chargesum(ineig(in,ia))
              end do
              if (abs(chargesum(ia)) .ge. qmin) then
                indexa(ia)=-1
                nsaltats=nsaltats+1
              end if
            end if
          end do
        end if
      end if
      if (ansrun .eq. 'l' .or. ansrun .eq. 'b') then
c       Atom selection by name
        lname=inamcol2-inamcol1+1
100     call getnamelist(saltatname,lname,nsaltatnams,
     -    'Salt-bridge forming atom names',30,100)
        print *,'NOTE: only oxygens will be assumend to be negative'
        do ia=1,nslt
          atnam=line(index(ia))(inamcol1:inamcol2)
          call leftadjustn(atnam,atnam,8)
          do in=1,nsaltatnams
            if (atnam(1:lname) .eq. saltatname(in)(1:lname)) then
              if (iatnum(ia) .ne. 1 .and. iatnum(ia) .ne. 6) then
                if (indexa(ia) .eq. 0) then
                  nsaltats=nsaltats+1
                  indexa(ia)=-1
                end if
              else
                write (6,2148) ia,atnam(1:lname),iatnum(ia)
                go to 100
              end if
            end if
c           Only oxygen is allowed to be negative
            if (iatnum(ia) .eq. 8) chargesum(ia)=-1.0
          end do
        end do
      end if
      if (nsaltats .eq. 0) then
        print *,'Solute contains no salt-bridge atoms'
        ifail=1
        return
      else
        write (6,2142) nsaltats
        write (iw0,2142) nsaltats
      end if
      nanchorr=0
      call zeroiti(indexn,0,nslt)
9134  ansrun=' '
      iall=0
      do while (iall .eq. 0 .and. ansrun .ne. 'q')
        call quiz(ansrun,iansrun,' ',' ',1,
     -    'hydrophobic/salt bridge/contact anchor atoms',44,0,5,6,75)
        ifail=0
        call definelist(ansrun,nslt,nanchorr,indexn,indexo,nsegslt,
     -    segid4,iresno,molsltlim,'salt bridge',11,iall,maxanchorlist)
      end do
      nanchor=0
c     indexa(ia) < 0: potential SB former
c     indexn(ja): atom selected for anchor
      do ja=1,nanchorr
        ia=indexn(ja)
        if (ia .gt. 0) then
          if (indexa(ia) .lt. 0) then
            indexa(ia)=1
            nanchor=nanchor+1
          end if
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
      write (iw0,*) 'Salt-bridge anchor atoms:'
      ia=0
      if (nanchor .gt. maxanchorlist) then
        print *,'Anchor list limit (',maxanchorlist,') is reached'
        ifail=1
      end if
      do ib=1,nslt
        if (indexa(ib) .gt. 0) then
          ia=ia+1
          ianchor(ia)=ib
          write (iw0,2146) ia,line(index(ib))(irescol1:irescol2),
     -      iresno(ib),segid4(isegno(ib)),
     -      line(index(ib))(inamcol1:inamcol2),ib,chargesum(ib)
        end if
        if (chargesum(ib) .lt. 0.0) then
          if (indexa(ib) .eq. 1) indexa(ib)=2
          if (indexa(ib) .eq. -1) indexa(ib)=-2
        end if
      end do
      if (ia .ne. nanchor) then
        print *,'PRORAM ERROR: ia=',ia,' nanchor=',nanchor
      end if
      return
2142  format(' Number of salt-bridge forming atoms=',i4)
2146  format(i5,' Residue=',a,' (',i5,') Chain/seg=',a,' Atom=',a,
     -  ' (',i6,') q=',f6.3)
2148  format(' Atomic number of atom ',i6,' (',a,')=',i3,/,
     -  ' Hydrogens and carbons are not allowed to form salt bridges')
      end
