      subroutine modify(inpcrdtyp,c,n,nslt,line,index,iatno,nneig,ineig,
     -  isegno,iresno,charge,icharges,natomadd,lineread,pi,maxneig,
     -  maxrec)
      dimension c(3,maxrec),index(maxrec),iatno(maxrec),nneig(maxrec),
     -  ineig(maxneig,maxrec),isegno(maxrec),iresno(maxrec),
     -  charge(maxrec)
      character*4 atnam,pflab
      character* 132 line(maxrec)
      character*2 iatnm2
      common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
     -  mmatno(64),iatnm2(99)
      character*5 crdext
      common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
     -  iommod,iommc,iommc4,iogro,iomol2,iomae,iocif,ioxxx,ioins,
     -  ionxyz,iosxyz,iosxyzrq,iograsp,iofull,lext(19),crdext(19)
      character*1 ans,modtyp
      character*2 ATD4typ
      character*80 question
      dimension c12(3),iadef(3)
      data ia4 /0/
      character*21 ngdefnum
      data ngdefnum /'Defining neighbour # '/,ATD4typ /'  '/
      if (inpcrdtyp .eq. iomol2) then
        print *,'Sorry, molecule editing with mol2 format is not ',
     -    'implemented'
         return
      end if
      call setcol(inpcrdtyp,ncol,idcol,ialtcol,iinscol,
     -  inamcol1,inamcol2,irescol1,irescol2,iccol1,iccol2,
     -  iresncol1,iresncol2,iseqncol1,iseqncol2,isegcol1,isegcol2,
     -  iresidcol1,iresidcol2,iqcol1,iqcol2,ipotcol1,ipotcol2,
     -  iocccol1,iocccol2,ichemcol1,ichemcol2,nrescol,nresncol,
     -  nsegcol,nnamcol,iofull)
      chargech=0.0
      iatnoadd=0
      do while (.true.)
        ifail=0
        atnam='    '
        call quiz(modtyp,iansmut,' ',' ',0,'modification type',17,
     -    0,5,6,89)
        if (modtyp .eq. 'm') then
          call getint('Atom number to replace',22,0,1,nslt,ia,0)
          call listatom(line,index,iatno,ia,inpcrdtyp,iofull,
     -      'Atom to mutate',14,n,maxrec)
          question='Name'//' of the new atom'
          lq=20
          call getname(atnam,len,question,lq,4,'',0,0,0,3)
          iatno(ia)=ianum(atnam,1,nnamcol)
          if (nneig(ia) .eq. 1) then
            ia1=ineig(1,ia)
            call listatom(line,index,iatno,ia1,inpcrdtyp,iofull,
     -        'Root atom     ',14,n,maxrec)
            call getint('Valence',7,nval(iatno(ia)),1,8,newn,0)
            call defaultbondl(newn,nneig(ia1),iatno(ia),iatno(ia1),
     -        r12)
            call getreal('New distance',12,r12,rij,1,0)
            call arrdiff(c(1,ia),c(1,ia1),c12,3)
            rnorm=sqrt(scprod(c12,c12))
            do k=1,3
              c(k,ia)=c(k,ia1)+rij*c12(k)/rnorm
            end do
          else
            write (6,1001) ia
          end if
        else if (modtyp .eq. 'a' .or. modtyp .eq. 'b' .or.
     -           modtyp .eq. 't') then
          call getint(
     -      'Index of the atom that the new atom is bonded to',48,
     -      0,1,nslt,ia1,0)
          call listatom(line,index,iatno,ia1,inpcrdtyp,iofull,
     -      'Atom R1         ',16,n,maxrec)
          question='Name'//' of the new atom'
          lq=20
          call getname(atnam,len,question,lq,4,'',0,0,0,3)
          iatnoadd=ianum(atnam,1,nnamcol)
          if (iatno(ia1) .eq. 1) then
            write (6,1002) ia1
            ifail=1
          end if
          if (modtyp .eq. 'a') then
            write (6,1000) atnam(1:len)
            call findroot('R2',ia1,0,iatno,nneig,ineig,line,index,
     -        inpcrdtyp,iofull,ia2,n,ifail,maxneig,maxrec)
            call findroot('R3',ia2,ia1,iatno,nneig,ineig,line,index,
     -        inpcrdtyp,iofull,ia3,n,ifail,maxneig,maxrec)
            if (ifail .eq. 0) then
              call getint('Valence',7,nval(iatnoadd),1,8,newn,0)
              call defaultbondl(nneig(ia1)+1,newn,iatno(ia1),iatnoadd,
     -          r12)
              call getreal('R1-X distance',13,r12,rij,1,0)
100           call quiz(ans,ians,' ',' ',0,'R2-R1-X angle',13,0,5,6,0)
              if (ans .eq. 'i' .or. ans .eq. 'I') then
                call getreal('R2-R1-X angle',13,999999.0,aijk,1,0)
              else if (ans .eq. '3') then
                aijk=180.0-acos(1.0/3.0)*180.0/pi
              else if (ans .eq. '2') then
                aijk=120.0
              else if (ans .eq. '1') then
                aijk=180.0
              else
                print *,'Invalid answer'
                go to 100
              end if
200           call quiz(ans,ians,' ',' ',0,
     -          'R3-R2-R1-X torsion angle',24,0,5,6,0)
              if (ans .eq. 'i' .or. ans .eq. 'I') then
                call getreal('R3-R2-R1-X torsion',18,999999.0,tijkl,0,0)
              else if (ans .eq. 'c' .or. ans .eq. 'C') then
                tijkl=0.0
              else if (ans .eq. 't' .or. ans .eq. 'T') then
                tijkl=180.0
              else if (ans .eq. '+') then
                tijkl=60.0
              else if (ans .eq. '-') then
                tijkl=-60.0
              else if (ans .eq. 'p' .or. ans .eq. 'P') then
                tijkl=120.0
              else if (ans .eq. 'm' .or. ans .eq. 'M') then
                tijkl=-120.0
              else
                print *,'Invalid answer'
                go to 200
              end if
              ia=n+1
              call addatom(1,c(1,ia3),c(1,ia3),c(1,ia2),c(1,ia1),
     -          c(1,ia),rij,aijk,tijkl,0.0,ia,pi,1,ifail)
            end if
          else if (modtyp .eq. 'b' .or. modtyp .eq. 't') then
            if (modtyp .eq. 'b') nndef=2
            if (modtyp .eq. 't') nndef=3
            call zeroiti(iadef,0,3)
            if (nneig(ia1) .ne. nndef)
     -        write (6,1007) ia1,nneig(ia1),nndef
            do ind=1,nndef
              write (ngdefnum(21:21),1003) ind
              if (ind .le. nneig(ia1)) then
c               Neighbor exists
                iadef(ind)=ineig(ind,ia1)
                call listatom(line,index,iatno,iadef(ind),inpcrdtyp,
     -            iofull,ngdefnum,21,n,maxrec)
               else
C               No bonded neighbor - ask for one
                write (6,1004) ind
                call getint('Index of the next definig atom',30,0,1,
     -            nslt,iadef(ind),0)
                call listatom(line,index,iatno,iadef(ind),inpcrdtyp,
     -            ioins,ngdefnum,21,n,maxrec)
               end if
            end do
            ia2=iadef(1)
            ia3=iadef(2)
            ia4=iadef(3)
            if (ia4 .eq. 0) ia4=ia3
            call getint('Valence',7,nval(iatnoadd),1,8,newn,0)
            call defaultbondl(nneig(ia1)+1,newn,iatno(ia1),iatnoadd,r12)
            call getreal('R1-X distance',13,r12,rij,1,0)
            bend=0.0
            if (nndef .eq. 2) call getreal(
     -        'Angle between R1-X and the negative bisector',44,0.0,
     -        bend,1,000)
            ia=n+1
            call addatom(nndef,c(1,ia4),c(1,ia3),c(1,ia2),c(1,ia1),
     -        c(1,ia),rij,aijk,tijkl,bend,ia,pi,1,ifail)
          end if
          if (ifail .eq. 0) then
            natomadd=natomadd+1
            if (lineread .gt. index(n)) then
c             Make room
              do im=lineread+1,index(n)+2,-1
                line(im)=line(im-1)
              end do
              lineread=lineread+1
            end if
            n=n+1
            nslt=nslt+1
            index(n)=index(n-1)+1
            line(index(n))=line(index(ia1))
            nneig(n)=1
            ineig(1,n)=ia1
            nneig(ia1)=nneig(ia1)+1
            ineig(nneig(ia1),ia1)=n
            isegno(n)=isegno(ia1)
            iresno(n)=iresno(ia1)
            call blankout(line(index(n)),iresncol1,iresncol2)
            icolw=iresncol1
            lenw=iresncol2-iresncol1+1
            call writeint(line(index(n)),icolw,iresno(n),lenw)
            call rightadjustline(line(index(n)),iresncol1,iresncol2)
            if (iseqncol2 .gt. iseqncol1) then
              call blankout(line(index(n)),iseqncol1,iseqncol2)
              icolw=iseqncol1
              lenw=iseqncol2-iseqncol1+1
              read (line(index(n-1))(iseqncol1:iseqncol2),*) lastseq
              call writeint(line(index(n)),icolw,lastseq+1,lenw)
              call rightadjustline(line(index(n)),iseqncol1,iseqncol2)
            end if
            if (line(index(n))(77:78) .ne. '  ' .and.
     -          iatnoadd .gt. 0) line(index(n))(77:78)=iatnm2(iatnoadd)
          end if
        else if (modtyp .eq. 'h') then
c         Add amide hydrogens

        else if (modtyp .eq. 'c') then
C         Create new bond
          call getint('Index of the first atom of the new bond',39,
     -      0,1,nslt,ia1,0)
          call listatom(line,index,iatno,ia1,inpcrdtyp,iofull,
     -          ' ',1,n,maxrec)
          call getint('Index of the second atom of the new bond',40,
     -      0,1,nslt,ia2,0)
          call listatom(line,index,iatno,ia2,inpcrdtyp,iofull,
     -          ' ',1,n,maxrec)
          ibonded=0
          do ia=1,nneig(ia1)
            if (ia2 .eq. ineig(ia,ia1)) ibonded=1
          end do
          if (ibonded .eq. 1) then
            print *,'These atoms are already bonded'
          else
            nneig(ia1)=nneig(ia1)+1
            nneig(ia2)=nneig(ia2)+1
            ineig(nneig(ia1),ia1)=ia2
            ineig(nneig(ia2),ia2)=ia1
          end if
        else if (modtyp .eq. 'q') then
          if (natomadd .gt. 0) then
            print *,natomadd,' atoms were added'
          end if
          if (icharges .gt. 0) write (6,1006) chargech
          return
        end if
        if (ifail .eq. 0) then
c         Specify properties of the new atom
          line(index(ia))(inamcol1:inamcol2)=atnam(1:nnamcol)
          if (ipotcol1 .le. ipotcol2) then
            question='Potential label'//' of the new atom'
            lq=31
            pflab='    '
            call getname(pflab,len,question,lq,ipotcol2-ipotcol1+1,'',0,
     -        0,0,0)
            line(index(ia))(ipotcol1:ipotcol2)=
     -        pflab(1:ipotcol2-ipotcol1+1)
          end if
          if (icharges .gt. 0) then
            if (modtyp .eq. 'm') then
              read (line(index(ia))(iqcol1:iqcol2),*) qprev
              chargech=chargech-qprev
            end if
            call getreal('Charge',6,999999.0,q,0,0)
            if (q .lt. 0.0) q=q-0.000001
            if (q .gt. 0.0) q=q+0.000001
            call putreal(line(index(ia))(iqcol1:iqcol2),iqcol2-iqcol1+1,
     -        q,4)
            chargech=chargech+q
            charge(ia)=q
          end if
          if (inpcrdtyp .eq.  ioa4pdb) then
            call getname(ATD4typ,len,'Autodock-4 atom type',20,2,'',0,0,
     -        0,0)
            line(index(ia))(78:79)=ATD4typ
          end if
          if (inpcrdtyp .eq. iocha) then
            icinc=0
            if (n .ge. 1000000) icinc=1
            write (line(index(ia))(iccol1+icinc:iccol2+icinc),1011)
     -        (c(k,ia),k=1,3)
          else if (inpcrdtyp .eq. iochaex) then
            write (line(index(ia))(iccol1:iccol2),1015) (c(k,ia),k=1,3)
          else if (inpcrdtyp .eq. iommc) then
            write (line(index(ia))(iccol1:iccol2),1011) (c(k,ia),k=1,3)
          else if (ispdb(inpcrdtyp) .gt. 0) then
            write (line(index(ia))(iccol1:iccol2),1012) (c(k,ia),k=1,3)
          else if (inpcrdtyp .eq. iommod) then
            write (line(index(ia))(iccol1:iccol2),1013) (c(k,ia),k=1,3)
          else if (inpcrdtyp .eq. iogro) then
            write (line(index(ia))(iccol1:iccol2),1012)
     -        (c(k,ia)/10.0,k=1,3)
          else if (inpcrdtyp .eq. ioins) then
            write (line(index(ia))(iccol1:iccol2),1014) (c(k,ia),k=1,3)
          end if
        end if
      end do
1000  format(' The position of the new atom ',a,' will be specified by',
     -  ' specifing ',/,' a sequence of bonded atoms R3-R2-R1-X,',
     -  ' the distance r(R1-X),',/,
     -  ' the angle a(R2-R1-X) and the torsion t(R3-R2-R1-X)')
1001  format(' Atom',i5,' has more than one bond - position will not ',
     -  'change')
1002  format(' Atom',i5,' is a hydrogen - can not attach to it')
1003  format(i1)
1004  format(' Neighbor #',i1,' is undefined')

1006  format(' Total charge changed by ',f7.4)
1007  format(' NOTE: atom',i6,' has',i2,' bonded neighbours ',
     -  '(instead of',i2,')')
1011  format(3f10.5)
1012  format(3f8.3)
1013  format(3f12.5)
1014  format(3f15.9)
1015  format(3f10.10)
      end
