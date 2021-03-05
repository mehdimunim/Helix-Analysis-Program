      subroutine rrdist(c,n,nres,ixresno,iresno,ixres,resnames,ifres,
     -  ilres,increst,incsolvrr,irefres1,irefres2,irefseg1,irefseg2,
     -  irefresinc,resdistlim,inegres1,inegres2,inegseg1,inegseg2,
     -  inegresinc,resapplim,listrefres,nrefres,nrefrange,listnegres,
     -  nnegres,nnegrange,iatnum,ignoreh,nframe,itypavg,iwrrdr,
     -  iwrrdc,iwrrmp,ir1,ir2,ic1,ic2,is1,is2,irescount1,irescount2,
     -  irescount3,line,index,indexn,indexo,itemp1,itemp2,itemp3,
     -  itemp4,temp1,temp2,maxrec,maxrsd)
      dimension c(3,n),iresno(n),ixresno(maxrsd),ixres(n),index(n),
     -  indexn(n),indexo(n),ifres(n),ilres(n),iatnum(n),
     -  listrefres(maxrsd),listnegres(maxrsd),
     -  irescount1(maxrsd),irescount2(maxrsd),irescount3(maxrsd),
     -  itemp1(n),itemp2(n),itemp3(n),itemp4(n),temp1(n),temp2(n)
      character*8 resnames(maxrsd)
      character* 132 line(maxrec)
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (6*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ sdm(MAX2D,MAX2D),nng(MAX2D),rmsd(MAX2D,MAX2D),
     -  ihbtores(MAXBONDS),nusepair(MAXBONDS),nhb_atot(MAXBONDS),
     -  nhb_rtot(MAXBONDS),a1(MAXBONDS),a2(MAXBONDS),fill(IFILL2)
      call header_rrdist(iwrrdr,nrefrange,nnegrange,resdistlim,
     -  irefseg1,ixresno(irefres1),irefseg2,ixresno(irefres2),irefres1,
     -  irefres2,inegseg1,ixresno(inegres1),inegseg2,ixresno(inegres2),
     -  inegres1,inegres2,incsolvrr,listrefres,listnegres,nrefres,
     -  nnegres,itemp1,maxrsd)
      call header_rrdist(iwrrdc,nrefrange,nnegrange,resapplim,
     -  irefseg1,ixresno(irefres1),irefseg2,ixresno(irefres2),irefres1,
     -  irefres2,inegseg1,ixresno(inegres1),inegseg2,ixresno(inegres2),
     -  inegres1,inegres2,incsolvrr,listrefres,listnegres,nrefres,
     -  nnegres,itemp1,maxrsd)
      call header_rrdist(iwrrmp,nrefrange,nnegrange,0.0,
     -  irefseg1,ixresno(irefres1),irefseg2,ixresno(irefres2),irefres1,
     -  irefres2,inegseg1,ixresno(inegres1),inegseg2,ixresno(inegres2),
     -  inegres1,inegres2,incsolvrr,listrefres,listnegres,nrefres,
     -  nnegres,itemp1,maxrsd)
      write (iwrrdr,1004) 'reference','neighbour','representative atoms'
      write (iwrrdc,1004) 'reference','neighbour','closest approach'
      if (ignoreh .eq. 1) write (iwrrdc,1002) 'Closest approach'
      if (nframe .eq. 1 .and. itypavg .gt. 0) then
        do irr=irefres1,irefres2
          do inr=inegres1,inegres2
            rmsd(irr-irefres1+1,inr-inegres1+1)=0.0
            sdm(irr-irefres1+1,inr-inegres1+1)=0.0
          end do
        end do
      end if
      resdlim2=resdistlim**2
      resalim2=resapplim**2
      call zeroiti(indexo,0,n)
      call zeroiti(indexn,0,n)
      call zeroiti(itemp1,0,n)
      call zeroiti(itemp2,0,n)
      irramin=0
      ncontact1=0
      distcontactsum1=0.0
      ncontact2=0
      distcontactsum2=0.0
      do irrr=1,nrefres
        irr=listrefres(irrr)
c       First find the representative atom and atom range for residue irr
        call findat(irarep,ifres(irr),ilres(irr),line,index,
     -    ir1,ic1,maxrec)
c       write (iwrrdc,7711) irr,irarep,irra1,irra2,
c    -     line(index(irarep))(ir1:ir1+3),
c    -     line(index(irarep))(ic1:ic2)
c7711  format('irr,irarep=',2i5,' irra1,irra2=',2i5,' r,a=',2a5)
        inramin=0
        do inrr=1,nnegres
          inr=listnegres(inrr)
          if (inr .ne. irr) then
            call findat(inarep,ifres(inr),ilres(inr),line,index,
     -        ir1,ic1,maxrec)
c           Calculate first based on representative atoms
            r2=dist2(c(1,irarep),c(1,inarep))

            if (r2 .lt. resdlim2) then
              indexo(irr)=irarep
              indexo(inr)=inarep
              call writeprox(iwrrdr,irarep,inarep,line,index,r2,
     -          iresno,ir1,ir2,ic1,ic2,is1,is2,maxrec)
              itemp1(irr)=1
              itemp1(inr)=1
              irescount1(irr)=irescount1(irr)+1
              irescount1(inr)=irescount1(inr)+1
              ncontact1=ncontact1+1
              distcontactsum1=distcontactsum1+r2
            end if
            if (itypavg .eq. 1) then
              rmsd(irr-irefres1+1,inr-inegres1+1)=
     -          rmsd(irr-irefres1+1,inr-inegres1+1)+sqrt(r2)
              sdm(irr-irefres1+1,inr-inegres1+1)=
     -          sdm(irr-irefres1+1,inr-inegres1+1)+r2
            end if
c           Now obtain closest approach
            call findapproach(c,ifres(irr),ilres(irr),ifres(inr),
     -        ilres(inr),iatnum,ignoreh,irarepm,inarepm,rmin2,maxrec)
            if (rmin2 .lt. resalim2) then
              indexn(irr)=1
              indexn(inr)=1
              call writeprox(iwrrdc,irarepm,inarepm,line,index,
     -          rmin2,iresno,ir1,ir2,ic1,ic2,is1,is2,maxrec)
              itemp2(irr)=1
              itemp2(inr)=1
              irescount2(irr)=irescount2(irr)+1
              irescount2(inr)=irescount2(inr)+1
              ncontact2=ncontact2+1
              distcontactsum2=distcontactsum2+rmin2
            end if
            if (itypavg .eq. 2) then
              rmsd(irr-irefres1+1,inr-inegres1+1)=
     -          rmsd(irr-irefres1+1,inr-inegres1+1)+sqrt(rmin2)
              sdm(irr-irefres1+1,inr-inegres1+1)=
     -          sdm(irr-irefres1+1,inr-inegres1+1)+rmin2
            end if
          end if
        end do
      end do
      if (ncontact1 .gt. 0)
     -  write (iwrrdr,1001) ncontact1,distcontactsum1/ncontact1
      if (ncontact2 .gt. 0)
     -  write (iwrrdc,1001) ncontact2,distcontactsum2/ncontact2
      call writeuniquelist(itemp1,ixresno,nres,resnames,ir2-ir1+1,
     -  iwrrdr,irefres1,irefres2,irefresinc,itemp3,itemp4,
     -  'reference',9,maxrsd)
      call writeuniquelist(itemp1,ixresno,nres,resnames,ir2-ir1+1,
     -  iwrrdr,inegres1,inegres2,inegresinc,itemp3,itemp4,
     -  'neighbour',9,maxrsd)
      call writeuniquelist(itemp2,ixresno,nres,resnames,ir2-ir1+1,
     -  iwrrdc,irefres1,irefres2,irefresinc,itemp3,itemp4,
     -  'reference',9,maxrsd)
      call writeuniquelist(itemp2,ixresno,nres,resnames,ir2-ir1+1,
     -  iwrrdc,inegres1,inegres2,inegresinc,itemp3,itemp4,
     -  'neighbour',9,maxrsd)
      call zeroiti(itemp1,0,n)
      do i=1,nrefres
        itemp1(listrefres(i))=1
      end do
      do i=1,nnegres
        itemp1(listnegres(i))=1
      end do
      call masktolist(itemp2,itemp1,nres,ncomplneg,0)
c     itemp2 is the list of residues not on the neighbor and reference ist
      if (ncomplneg .gt. 0 .and. increst .gt. 0) then
        write (iwrrdr,1004) 'reference','other','representative atoms'
        write (iwrrdc,1004) 'reference','other','closest approach'
        irramin=0
        do irrr=1,nrefres
          irr=listrefres(irrr)
          if (indexo(irr) .gt. 0 .or. indexn(irr) .eq. 0) then
            call findat(irarep,ifres(irr),ilres(irr),line,index,
     -        ir1,ic1,maxrec)
            inramin=0
            do irrr1=1,ncomplneg
              irr1=itemp2(irrr1)
              if ((irr1 .lt. inegres1 .or. irr1 .gt. inegres2)
     -           .and. irr .ne. irr1) then
                call findat(inarep,ifres(irr1),ilres(irr1),line,
     -            index,ir1,ic1,maxrec)
                if (indexo(irr) .gt. 0 .and. indexo(irr1) .eq. 0) then
                  r2=dist2(c(1,irarep),c(1,inarep))
                  if (r2 .lt. resdlim2)
     -              call writeprox(iwrrdr,irarep,inarep,line,
     -                index,r2,iresno,ir1,ir2,ic1,ic2,is1,is2,maxrec)
                end if
                if (indexn(irr) .gt. 0 .and. indexn(irr1) .eq. 0) then
                  call findapproach(c,ifres(irr),ilres(irr),
     -              ifres(irr1),ilres(irr1),iatnum,ignoreh,irarepm,
     -              inarepm,r2,maxrec)
                  if (r2 .lt. resalim2)
     -              call writeprox(iwrrdc,irarepm,inarepm,line,
     -                index,r2,iresno,ir1,ir2,ic1,ic2,is1,is2,maxrec)
                  end if
                end if
            end do
          end if
        end do
        write (iwrrdr,1004) 'neighbour','other','representative atoms'
        write (iwrrdc,1004) 'neighbour','other','closest approach'
        inramin=0
        do irrr=1,nnegres
          irr=listnegres(irrr)
          if (indexo(irr) .gt. 0 .or. indexn(irr) .eq. 0) then
            call findat(inarep,ifres(irr),ilres(irr),line,index,
     -        ir1,ic1,maxrec)
            irramin=0
            do irrr1=1,ncomplneg
              irr1=itemp2(irrr1)
              if ((irr1 .lt. irefres1 .or. irr1 .gt. irefres2)
     -           .and. irr .ne. irr1) then
                call findat(irarep,ifres(irr1),ilres(irr1),line,
     -            index,ir1,ic1,maxrec)
                if (indexo(irr) .gt. 0 .and. indexo(irr1) .eq. 0) then
                  r2=dist2(c(1,irarep),c(1,inarep))
                  if (r2 .lt. resdlim2)
     -              call writeprox(iwrrdr,inarep,irarep,line,
     -                index,r2,iresno,ir1,ir2,ic1,ic2,is1,is2,maxrec)
                end if
                if (indexn(irr) .gt. 0 .and. indexn(irr1) .eq. 0) then
                  call findapproach(c,ifres(irr),ilres(irr),
     -              ifres(irr1),ilres(irr1),iatnum,ignoreh,irarepm,
     -              inarepm,r2,maxrec)
                  if (r2 .lt. resalim2)
     -              call writeprox(iwrrdc,irarepm,inarepm,line,
     -                index,r2,iresno,ir1,ir2,ic1,ic2,is1,is2,maxrec)
                  end if
                end if
            end do
          end if
        end do
      end if
c     Calculate contact pairs
      write (iwrrmp,1003)
      if (ignoreh .eq. 1) write (iwrrmp,1002) 'Contact calculation'
c     Put the ref and neg atom lists into itemp3, itemp4. resp
      nrefat=0
      do irr=1,nrefres
        ir=listrefres(irr)
        do ia=ifres(ir),ilres(ir)
          if (ignoreh .eq. 0 .or. iatnum(ia) .gt. 1) then
            nrefat=nrefat+1
            itemp3(nrefat)=ia
          end if
        end do
      end do
      nnegat=0
      do irr=1,nnegres
        ir=listnegres(irr)
        do ia=ifres(ir),ilres(ir)
          if (ignoreh .eq. 0 .or. iatnum(ia) .gt. 1) then
            nnegat=nnegat+1
            itemp4(nnegat)=ia
          end if
        end do
      end do
c     Find the proximal atoms for the ref and neg residues
      do ia=1,n
        temp1(ia)=100000.0
        temp2(ia)=100000.0
      end do
      do iaa=1,nrefat
        ia=itemp3(iaa)
        do jaa=1,nnegat
          ja=itemp4(jaa)
          rijsq=dist2(c(1,ia),c(1,ja))
          if (rijsq .lt. temp1(ia)) then
            itemp1(ia)=ja
            temp1(ia)=rijsq
          end if
          if (rijsq .lt. temp2(ja)) then
            itemp2(ja)=ia
            temp2(ja)=rijsq
          end if
        end do
      end do
      ncontact1=0
      distcontactsum1=0.0
      call zeroiti(indexo,0,n)
      call zeroiti(indexn,0,n)
      do iaa=1,nrefat
        ia=itemp3(iaa)
        if (itemp1(ia) .gt. 0) then
          if (itemp2(itemp1(ia)) .eq. ia) then
c           Mutually proximal atom pair found
            call writeprox(iwrrmp,ia,itemp1(ia),line,index,temp1(ia),
     -        iresno,ir1,ir2,ic1,ic2,is1,is2,maxrec)
            indexo(ixres(ia))=1
            indexn(ixres(itemp1(ia)))=1
            irescount3(ixres(ia))=irescount3(ixres(ia))+1
            irescount3(ixres(itemp1(ia)))=irescount3(ixres(ia))+1
              ncontact1=ncontact1+1
              distcontactsum1=distcontactsum1+temp1(ia)
          end if
        end if
      end do
      if (ncontact1 .gt. 0)
     -  write (iwrrmp,1001) ncontact1,distcontactsum1/ncontact1
      call writeuniquelist(indexo,ixresno,nres,resnames,ir2-ir1+1,
     -  iwrrmp,irefres1,irefres2,irefresinc,itemp1,itemp2,
     -  'reference',9,maxrsd)
      call writeuniquelist(indexn,ixresno,nres,resnames,ir2-ir1+1,
     -  iwrrmp,inegres1,inegres2,inegresinc,itemp1,itemp2,
     -  'neighbour',9,maxrsd)
      return
1001  format(' Number of residue pairs within threshold=',i6,
     -  ' Average distance=',f5.2,' A')
1002  format(1x,a,' is based on heavy atoms only')
1003  format(' Contact atom pairs are the atoms that are mutually ',
     -  'proximal to each other')
1004  format(/,' Residue distances between the ',a,' and ',a,
     -  ' residues, based on the',/,1x,a,':')
      end
