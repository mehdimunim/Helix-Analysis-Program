      subroutine writeconf(iout,inpcrdtyp,iotyp,inpcrdtyporg,n0,n,nslt,
     -  naslv,islvw,iasv,namesv,qsv,pflsv,icreaterec,iwhead,imodel,
     -  iatnum,ifchrg,nconfig,innlist,c,rprox,cv,ixres,iresno,atnames,
     -  resnames,segnames,charge,isegno,altcol,inscol,ninsres,marker,
     -  ntitlin,ntitlinw,title,ireseq,iresnrestart,iresidrestart,nneig,
     -  nneiga,nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,npneig,numres,
     -  numslv,resnamslv,line,blankline,mmtype,ibnd,index,indexn,indexo,
     -  iclean,molresflag,idupl,ihetat,hblimfac,angmin,noleftad,nosort,
     -  noreg,iusecvforq,nobondord,nqdec,iqspaceask,ianaltyp,ifromtraj,
     -  etot,ietot,nclstmem,noend,keeprem,iwriteatsym,radtodeg,
     -  maxrepconf,maxng,maxrsd,maxrec)
      character* 132 line(maxrec),blankline
      character*80 title
      character*1 altcol(maxrec),inscol(maxrec)
      character*4 segnames(maxrsd),pflsv(100)
      character*8 atnames(maxrec),resnames(maxrsd),namesv(100)
      character*6 marker(16)
      dimension nneig(maxrec),ineig(maxng,maxrec),iatnum(maxrec),
     -  ifchrg(maxrec),c(3,maxrec),cv(maxrec),ixres(maxrec),
     -  rprox(maxrec),charge(maxrec),nhbneig(maxrec),nneiga(maxrec),
     -  nhneig(maxrec),nnneig(maxrec),ncneig(maxrec),nsneig(maxrec),
     -  npneig(maxrec),isegno(maxrec),iresno(maxrec),mmtype(maxrec),
     -  ibnd(maxng,maxrec),index(maxrec),indexn(maxrec),indexo(maxrec),
     -  molresflag(maxrsd),idupl(maxrsd),ihetat(maxrec),iasv(100),
     -  qsv(100)
      common /columnlim/ incol(19),iidcol(19),iialtcol(19),iiinscol(19),
     -  iinamcol(2,19),iirescol(2,19),iiccol(2,19),iiresncol(2,19),
     -  iiseqncol(2,19),iisegcol(2,19),iiresidcol(2,19),iiqcol(2,19),
     -  iipotcol(2,19),iiocccol(2,19),iichemcol(2,19)
      character*11 formatname
      common /formats/ iqdconv(20),formatname(19)
      character*5 crdext
      common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
     -  iommod,iommc,iommc4,iogro,iomol2,iomae,iocif,ioxxx,ioins,
     -  ionxyz,iosxyz,iosxyzrq,iograsp,iofull,lext(19),crdext(19)
      character*2 iatnm2
      common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),mmatno(64),
     -  iatnm2(99)
      character*1 abc,digits,hexdigits
      common /charactersets/ ihex(25),abc(62),digits(14),hexdigits(25)
      common /logging/ logfile,ipredict
      character*2 mmcgm
      character*4 segnam,chnam,xseg
      character*6 potnam
      character*8 resnam,resnamslv,atomnam
      character*8 segnamlong
      character*31 question
      dimension iabc(62),cw(3)
      data xseg /'XXXX'/
c     print *,'WRITECONF nslt=',nslt,' n0,n=',n0,n,
c    -  ' iwriteatsym=',iwriteatsym
      mmcgm='  '
      call set_pdbinfo(iotyp,iwriteatsym,iwritecon,iobpdb,iocpdb,0)
      call setcol(inpcrdtyp,ncol,idcol,ialtcol,iinscol,
     -  inamcol1,inamcol2,irescol1,irescol2,iccol1,iccol2,
     -  iresncol1,iresncol2,iseqncol1,iseqncol2,isegcol1,isegcol2,
     -  iresidcol1,iresidcol2,iqcol1,iqcol2,ipotcol1,ipotcol2,
     -  iocccol1,iocccol2,ichemcol1,ichemcol2,nrescol,nresncol,
     -  nsegcol,nnamcol,iofull)
      frocc=1.0
      if (iotyp .eq. iommod .and. nobondord .eq. 0) then
c       Macromodel
        call nnlist(nslt,islvw,naslv,n,iatnum,ifchrg,c,nneig,nneiga,
     -    nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,npneig,line,
     -    irescol1,irescol2,inamcol1,inamcol2,index,nconfig,innlist,
     -    molresflag,hblimfac,angmin,0,ibnd,indexo,isegno,ixres,
     -    maxrepconf,0,0,radtodeg,0,maxng,maxng,maxrsd,maxrec)
        call bondord(iatnum,mmtype,n,nneig,ineig,nhneig,ibnd,maxng,c,
     -    index,ncneig,nsneig,inamcol1,inamcol2,irescol1,irescol2,line,
     -    nconfig,maxrepconf,maxrec)
      end if
      if (inpcrdtyp .le. ioins .and. isegcol2 .ge. isegcol1) then
c       Get segment id name list
        nsegm=isegno(n)
        segnames(1)='    '
        call getseg4(segnames(1),line(index(1)),isegcol1,nsegcol)
c       segnames(1)(1:nsegcol)=line(index(1))(isegcol1:isegcol2)
        nsegcolin=isegcol2-isegcol1+1
        nsegcolout=iisegcol(2,iotyp)-iisegcol(1,iotyp)+1
        isegnoprev=0
        do ia=1,min0(n0,n)
c         write(77,*) 'nconfig=',nconfig,' ia=',ia,'isegno=',isegno(ia)
c         write(77,*) 'index=',index(ia)
          if (isegno(ia) .ne. isegnoprev) then
            segnam='    '
            call getseg4(segnam,line(index(ia)),isegcol1,nsegcol)
c           segnam(1:nsegcol)=line(index(ia))(isegcol1:isegcol2)
            call leftadjust4(segnam,segnam)
            if (nsegcolout .lt. nsegcolin) then
              call lastchar(segnam,lc,4)
              if (lc .gt. 1 .or. nsegcolout .gt. 1) then
                question='New name of segment #'
                write (question(22:24),1005) isegno(ia)
                lq=27+lc
                question(25:lq)=' ('//segnam(1:lc)//')'
                lc=min0(lc,nsegcolout)
                if (lc .eq. 1) segnam(1:1)=abc(isegno(ia))
                call getname(segnam,len,question,lq,4,segnam(1:lc),lc,0,
     -            0,0)
              end if
            end if
            segnames(isegno(ia)-isegno(1)+1)=segnam
            isegnoprev=isegno(ia)
          end if
        end do
        if (n .gt. n0) then
          incsg=0
          if (n0 .eq. nslt) then
            segnames(isegno(n0)-isegno(1)+1)=xseg(1:min0(4,nsegcol))
            incsg=1
          end if
          do ia=n0+1,n
            isegno(ia)=isegno(n0)+incsg
          end do
        end if
        if (ischarmm(iotyp) .eq. 1 .or. iotyp .eq. iocpdb) then
          do is=1,nsegm
            segnam=segnames(is)
            if (nconfig .gt. 1) then
c             For Charmm output, change the segment ID when there is room
CHX
              do ic=1,4
                if (segnam(ic:ic) .eq. ' ') then
                  nslen=ic-1
                  go to 100
                end if
              end do
c             No blank, leave it
              go to 101
100           write (segnamlong,1008) nconfig
              call leftadjustn(segnamlong,segnamlong,8)
              call nextchar(segnamlong,ic,8)
              nnlen=min0(4-nslen,8-ic+1)
              segnam(nslen+1:nslen+nnlen)=segnamlong(ic:ic+nnlen-1)
101           continue
            else if (nsegcol .eq. 1) then
              segnam(2:4)='CHR'
            end if
            segnames(is)=segnam
          end do
        else if (iotyp .eq. iobpdb .and. nsegcol .gt. 1) then
c         Make sure new segment ID's are different
          call zeroiti(iabc,0,62)
          ndupl=0
          do is=1,nsegm
            do iss=1,62
              if (segnames(is)(1:1) .eq. abc(iss)) then
                if (iabc(iss) .eq. 0) then
                  iabc(iss)=is
                else
                  ndupl=ndupl+1
                  idupl(ndupl)=is
                end if
                go to 200
              end if
            end do
            if (segnames(is)(1:1) .eq. ' ') then
              ndupl=ndupl+1
              idupl(ndupl)=is
            else
              print *,'Invalid chain id character:',segnames(is)(1:1)
              segnames(is)(1:1)=' '
            end if
200         continue
          end do
          if (ndupl .gt. 0) then
            ifree0=0
            do id=1,ndupl
              do is=ifree0+1,62
                if (iabc(is) .eq. 0) then
                  ifree0=is
                  segnames(idupl(id))=abc(is)
                  go to 201
                end if
              end do
            end do
            ifree0=62
            print *,'More than 62 different segments exist - all the ',
     -        'remainings will be called Z'
201         continue
          end if
        end if
      else
c       Create segment id's from segment numbers
        do is=isegno(1),isegno(n)
          segnames(is)(1:1)=abc(min0(62,is))
          segnames(is)(2:4)='   '
        end do
      end if
c     Create atom records
      segnam='A   '
      resnam='RES     '
      atomnam='        '
      potnam='      '
      chnam='    '
      q=0.0
      icinc=0
      call zeroiti(ihetat,0,n)
      do iat=1,n
        if (line(index(iat))(1:6) .eq. 'HETATM') ihetat(iat)=1
      end do
      do iat=1,n
        if (inpcrdtyp .le. ioins .and. iat .le. n0) then
c         Full information is available
          if (isegcol2 .ge. isegcol1)
     -      segnam=segnames(isegno(iat)-isegno(1)+1)
          resnam(1:nrescol)=resnames(ixres(iat))(1:nrescol)
c         if (iat .lt. 25) print *,'RESNAM=',resnam
c    -      line(index(iat))(irescol1+icinc:irescol2+icinc)
          call leftadjustn(resnam,resnam,8)
          ic1=iresncol1
          call nextchar(line(index(iat)),ic1,132)
          ic2=ic1
          call nextblank(line(index(iat)),ic2,132)
          ic2=min0(iresncol2+1,ic2)
          if (ic2 .lt. ic1) then
            write (6,1000) iat
            ires=0
          else
c           call readint(line(index(iat)),ic1,ic2-1,ires,2,1,irerr)
            ires=iresno(iat)
          end if
          atomnam(1:nnamcol)=atnames(iat)(1:nnamcol)
          if (iotyp .eq. iommod) then
               write (potnam(1:4),1006) mmtype(iat)
          else if (ipotcol2 .ge. ipotcol1) then
             potnam(1:ipotcol2-ipotcol1+1)=
     -         line(index(iat))(ipotcol1:ipotcol2)
          end if
          if (iqcol2 .ge. iqcol1) q=charge(iat)
          if (iat .eq. 999999 .and. inpcrdtyp .eq. iocha) icinc=1
        else if (inpcrdtyp .le. ioins .and. iat .gt. n0) then
c         Add solvent information
          resnam=resnamslv
          iatslv=mod(iat-n0-1,naslv)+1
          iatnum(iat)=iasv(iatslv)
          charge(iat)=qsv(iatslv)
          q=qsv(iatslv)
          if (inpcrdtyp .eq. ioins .or. inpcrdtyp .eq. iommc)
     -      potnam(1:4)=pflsv(iatslv)
          atomnam=namesv(iatslv)
          ires=iresno(n0)+(iat-n0-1)/naslv+1
          iresno(iat)=ires
        else
c         Only limited information is available
          atomnam(1:2)=iatnm2(iatnum(iat))
          if (inpcrdtyp .eq. iosxyzrq) then
            ires=iresno(iat)
            q=charge(iat)
          end if
        end if
        if (iusecvforq .eq. 1) then
          q=cv(iat)
          frocc=rprox(iat)
        end if
        if (inpcrdtyp .eq. iocif) then
          read (line(index(iat))(iocccol1:iocccol2),*) frocc
          read (line(index(iat))(iqcol1:iqcol2),*) q
        end if
        chnam(1:2)=iatnm2(iatnum(iat))
        if (icreaterec .eq. 1) then
          irec=ntitlin+iat
          if (inpcrdtyp .eq. iotyp .and. iotyp .eq. iommc)
     -      mmcgm=line(index(iat))(51:52)
          line(irec)=blankline
          iqspace=1
          if (q .ge. 1000.0) then
            if (iqspaceask .eq. 1) then
              print *,'B-factor column entry > 1000 found'
              if (ipredict .eq. 0) then
                call askyn('Do you want to use all 6 characters',35,0,
     -            -1,iqspace,124,0)
              else
                print *,'All 6 characters will be used'
                iqspace=0
              end if
              iqspaceask=0
            end if
          end if
c          write (77,7611) iat,ires,q,segnam
c7611      format(' iat=',i6,' ires=',i5,' q=',f10.5,' segnam=',a)
          if (iotyp .le. ioins) then
            call trnsfr(cw,c(1,iat),3)
            if (inpcrdtyp .eq. iogro .and. iotyp .ne. iogro) then
              do k=1,3
                cw(k)=c(k,iat)*10.0
              end do
            else if (iotyp .eq. iogro) then
              do k=1,3
                cw(k)=c(k,iat)*0.1
              end do
            end if
            if (inpcrdtyp .eq. iomae) then
              segnam(1:1)=abc(min0(62,isegno(iat)))
              segnam(2:4)='   '
            end if
            call createrec(line(irec),inpcrdtyp,iotyp,cw(1),cw(2),cw(3),
     -        altcol(iat),inscol(iat),atomnam,resnam,segnam,iat,ires,
     -        ires,chnam,potnam,frocc,q,nqdec,iqspace,mmcgm,nneig(iat),
     -        ineig(1,iat),nconfig,ibnd(1,iat),ihetat(iat),blankline)
          else
            call writefree(line(irec),iotyp,chnam,c(1,iat),c(2,iat),
     -        c(3,iat),ires,q)
          end if
          index(iat)=irec
        end if
c       Left justify residue and atom names
        if (iotyp .eq. iogro) then
          call leftadjustline(line(index(iat)),iirescol(1,iotyp),
     -      iirescol(2,iotyp))
          call rightadjustline(line(index(iat)),
     -      iinamcol(1,iotyp),iinamcol(2,iotyp))
        else
          call leftadjustline(line(index(iat)),iirescol(1,iotyp),
     -      iirescol(2,iotyp))
          if ((ispdb(iotyp) .eq. 0 .or. inpcrdtyp .eq. iogro) .and.
     -         noleftad .eq. 0)
     -      call leftadjustline(line(index(iat)),
     -        iinamcol(1,iotyp),iinamcol(2,iotyp))
        end if
        if (iotyp .eq. ioins) call leftadjustline(line(index(iat)),
     -    iiresncol(1,iotyp),iiresncol(2,iotyp))
        if (ispdb(iotyp) .gt. 0 .and. noreg .eq. 0) call regularpdb(
     -    line(index(iat))(iinamcol(1,iotyp):iinamcol(2,iotyp)),
     -    line(index(iat))(iinamcol(1,iotyp):iinamcol(2,iotyp)),1)
      end do
      if (iclean .eq. 1) then
        if (iotyp .lt. ioins .and. nosort .eq. 0)
     -    call sortatres(line,n,index,indexn,indexo,nhbneig,nneiga,
     -    ncneig,nsneig,npneig,iiresncol(1,iotyp),iiresncol(2,iotyp),
     -    isegno,'      ',iotyp,nconfig,iout,maxrepconf,maxrec)
        if ((ischarmm(iotyp) .eq. 1 .or. ispdb(iotyp) .gt. 0 .or.
     -       iotyp .eq. iommod .or. iotyp .eq. iommc) .and.
     -       nosort .eq. 0)
     -    call reseq(line,index,n,nslt,nsegm,iatfirst,iresfirst,
     -      iresidfirst,iotyp,iresno,resnames,numres,numslv,naslv,
     -      ninsres,resnamslv,iresnrestart,iresidrestart,nconfig,ireseq,
     -      +1,maxrsd,maxrec)
      end if
      imodelnum=imodel
      if (imodel .gt. 0) imodelnum=nconfig
      call writeout(iout,inpcrdtyporg,iotyp,line,index,isegno,n,marker,
     -  iwhead,imodelnum,ntitlin,ntitlinw,title,blankline,nosort,
     -  ianaltyp,ifromtraj,etot,ietot,nclstmem,noend,keeprem,
     -  iwriteatsym,iatnum,maxrec)
      return
1000  format(' Atom ',i6,' has no residue number - it is set to zero')
1005  format(i3)
1006  format(i4)
1008  format(i8)
      end
