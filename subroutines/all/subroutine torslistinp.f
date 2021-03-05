      subroutine torslistinp(ixtor1234,talab,ltalab,ntang,inpcrdtyp,
     -  ioins,line,index,n,nslt,nresslt,iatnum,nneig,ineig,nhneig,ixres,
     -  ixresno,nsegm,isegno,ifres,ilres,iresno,indexs,listrefres,
     -  irescol1,irescol2,inamcol1,inamcol2,maxng,maxrsd,maxrec,maxtors)
      dimension ixtor1234(4,maxtors),ltalab(maxtors)
      dimension nneig(n),ineig(maxng,n),iatnum(n),nhneig(n),index(n),
     -  ixres(maxrec),ixresno(maxrsd),isegno(maxrec),ifres(maxrsd),
     -  ilres(maxrsd),iresno(maxrec),indexs(maxrec),listrefres(maxrsd)
      character*30 talab(maxtors)
      character* 132 line(maxrec)
      character*1 ansrun
      character*4 atnami,atnamj
      character*8 resnam
      character*38 question
      call quiz(ansrun,iansrun,' ',' ',1,'torsion list input',18,0,
     -    5,6,0)
      nohrot=0
      if (iansrun .gt. 3) call askyn(
     -    'Do you want to include torsions moving hydrogens only',53,
     -    0,-1,nohrot,0,0)
      ic1=inamcol1
      ic2=inamcol2
      ir1=irescol1
      ir2=irescol2
      nrescol=ir2-ir1+1
      call getresrange(nsegm,indexs,isegno,ixres,iresno,ifres,
     -  'residue to plot',15,nresslt,nslt,irefres1,irefres2,1,
     -  irefseg1,irefseg2,listrefres,nrefres,nrefrange,0,0,maxrsd,
     -  maxrec,111)
      if (iansrun .lt. 5) then
        ntang=0
c       print *,'RESNO RANGE:',irefres1,irefres2
c       print *,'ATOM RANGE:',ifres(irefres1),ilres(irefres2)
        do ia=ifres(irefres1),ilres(irefres2)
          resnam(1:nrescol)=line(index(ia))(ir1:ir2)
          call leftadjust4(line(index(ia))(ic1:ic2),atnami)
          do j=1,nneig(ia)
            ja=ineig(j,ia)
            nhneig0=nhneig(ja)*nohrot
            if (ia .lt. ja .and. nneig(ia)-nhneig(ia) .gt. 1 .and.
     -          nneig(ja)-nhneig0 .gt. 1) then
              call leftadjust4(line(index(ja))(ic1:ic2),atnamj)
c             Bond ia-ja found (on a heavy atom chain when nohrot=1)
              call checktorbond(resnam,ixres(ia),ixres(ja),atnami,
     -          atnamj,fixbond,ipep,ican,icac,issb)
              iskiptor=0
              if (ican .eq. 1 .and. ixres(ia) .eq. irefres1 .or.
     -            icac .eq. 1 .and. ixres(ia) .eq. irefres2) iskiptor=1
              ihfound=0
c             print *,'IA,ICAN,ICAC,IXRES=',ia,ican,icac,ixres(ia)
              if (iatnum(ia) .eq. 1 .or. iatnum(ja) .eq. 1) then
                write (6,2002) ia,resnam(1:nrescol),atnami,ja,
     -            line(index(ja))(ir1:ir2),atnamj
                ihfound=1
              end if
c             Loop and normal torsion stepsizes
              if (fixbond+ihfound+iskiptor .eq. 0) then
                dloop=30.0
                dtor=30.0
                if (iansrun .eq. 4
     -            .or. (iansrun .eq. 1 .and. ipep .eq. 0)
     -            .or. (iansrun .eq. 2 .and. ipep+ican+icac+issb .eq. 0)
     -            .or. (iansrun .eq. 3 .and. ican+icac .eq. 1)) then
                  if (ntang .eq. maxtors) then
                    write (6,1000) maxtors
                    go to 100
                  end if
                  ntang=ntang+1
c                 Bond ia-ja
                  do ii=1,nneig(ia)
                    iaa=ineig(ii,ia)
                    call leftadjust4(line(index(iaa))(ic1:ic2),atnami)
                    if (iatnum(iaa)*nohrot .ne. 1 .and.
     -                 (iansrun .ne. 3 .or. atnami .ne. 'CB  ')) then
                      ixtor1234(2,ntang)=ia
                      ixtor1234(3,ntang)=ja
                      if (ican .eq. 1) then
                        do jj=1,nneig(ia)
                          iaa=ineig(jj,ia)
                          if (iatnum(iaa) .eq. 6 .and. iaa .ne. ja)
     -                      ixtor1234(1,ntang)=iaa
                        end do
                        do jj=1,nneig(ja)
                          jaa=ineig(jj,ja)
                          call leftadjust4(line(index(jaa))(ic1:ic2),
     -                      atnamj)
                          if (atnamj(1:4) .eq. 'C   ')
     -                      ixtor1234(4,ntang)=jaa
                        end do
                      else if (icac .eq. 1) then
                        do jj=1,nneig(ia)
                          iaa=ineig(jj,ia)
                          if (iatnum(iaa) .eq. 7) ixtor1234(1,ntang)=iaa
                        end do
                        do jj=1,nneig(ja)
                          jaa=ineig(jj,ja)
                          if (iatnum(jaa) .eq. 7) ixtor1234(4,ntang)=jaa
                        end do
                      else
                        do jj=1,nneig(ja)
                          jaa=ineig(jj,ja)
                          if (iatnum(jaa)*nohrot .ne. 1 .and.
     -                        jaa .ne. ia .and. iaa .ne. ja .and.
     -                        iaa .ne. jaa) then
                            ixtor1234(1,ntang)=iaa
                            ixtor1234(4,ntang)=jaa
                          end if
                        end do
                      end if
                      if (ixtor1234(1,ntang)*ixtor1234(1,ntang) .eq. 0)
     -                  then 
                        print *,' Peptide bond C- or N+ is missing'
                        go to 100
                      end if
                    end if
                  end do
                end if
              end if
            end if
          end do
        end do
100      print *,'Number of torsion angles generated=',ntang
      else
        call getint('Number of torsions to track',27,0,1,49,ntang,0)
        do it=1,ntang
8004      question='Angle   , atomindices ( 4 numbers )'
          nerr=0
          write (question(6:8),2043) it
          call getintline(question,35,1,nslt,ixtor1234(1,it),4,0)
          if (inpcrdtyp .le. ioins)
     -      write (6,2039) 'T',(ixtor1234(k,it),
     -        line(index(ixtor1234(k,it)))(inamcol1:inamcol2),
     -        line(index(ixtor1234(k,it)))(irescol1:irescol2),k=1,4)
c         print *,' ixtor:', (ixtor1234(k,it),k=1,4)
          nerr=0
          do k=2,4
            if (isbonded(ixtor1234(k-1,it),ixtor1234(k,it),
     -          nneig,ineig,n,maxng) .eq. 0) then
              write (6,2066) (ixtor1234(kk,it),kk=k-1,k),' NOT '
              nerr=nerr+1
            end if
          end do
          do k=3,4
            if (isbonded(ixtor1234(k-2,it),ixtor1234(k,it),
     -          nneig,ineig,n,maxng) .eq. 1) then
              write (6,2066) (ixtor1234(kk,it),kk=k-1,k),' '
              nerr=nerr+1
            end if
          end do
          if (isbonded(ixtor1234(1,it),ixtor1234(4,it),
     -        nneig,ineig,n,maxng) .eq. 1) then
            write (6,2066) ixtor1234(1,it),ixtor1234(4,it),' '
            nerr=nerr+1
          end if
          if (nerr .gt. 0) then
            call askyn('Do you want to use this angle',29,1,-1,iok,
     -        0,0)
            if (iok. eq. 0) go to 8004
          end if
        end do
      end if
      do it=1,ntang
        call blankout(talab(it),1,30)
        if (inpcrdtyp .le. ioins) write (talab(it),2069)
     -    (line(index(ixtor1234(k,it)))(inamcol1:inamcol1+3),k=1,3),
     -    line(index(ixtor1234(3,it)))(irescol1:irescol2),
     -    ixresno(ixres(ixtor1234(3,it))),
     -    line(index(ixtor1234(4,it)))(inamcol1:inamcol1+3)
        ltalab(it)=30
      end do
      return
1000  format(' ERROR: maximum number of torsions (',i4,') exceeded',/,
     -  'Redimension Simulaid with increased value of the parameter ',
     -  'MAXCOPY')
2002  format(' WARNING: bond ',i5,' (',a,1x,a,') - ',i5,' (',a,1x,a,
     -  ') inlcudes a hydrogen',/,
     -  10x,'Torsion dropped - check the solute geometry')
2043  format(i3)
2039  format(1x,a1,':',i5,' (',a,1x,a,')',3(' -',i5,' (',a,1x,a,')'))
2066  format(' WARNING: atoms',i6,' and',i6,' are',a,'bonded',/,
     -  ' You can create bonds with the Edit option')
2069  format(2(a4,'-'),a4,'(',a4,i5,')-',a4)
      end
