      subroutine torsinp(nconfig,inpcrdtyp,itform,extnam,analfile,
     -  inpfile,namleni,line,index,n,nslt,c,iatnum,ifchrg,nneig,nneiga,
     -  nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,npneig,isegno,ixres,
     -  ibnd,indexo,naslv,islvw,innlist,molresflag,hblimfac,angmin,
     -  irescol1,irescol2,iresncol1,iresncol2,inamcol1,inamcol2,mask,
     -  itemp1,itemp2,pi,radtodeg,maxrepconf,maxng,maxbox,maxrsd,maxrec)
      dimension nneig(n),ineig(maxng,n),iatnum(n),ifchrg(n),c(3,n),
     -  nhbneig(n),nneiga(n),nhneig(n),nnneig(n),ncneig(n),nsneig(n),
     -  npneig(n),index(n),indexo(n),isegno(n),ixres(n),
     -  molresflag(maxrsd),ibnd(maxbox,maxrec),mask(maxrec),
     -  itemp1(maxrec),itemp2(maxrec)
      character*1 ansrun
      character*4 extnam,lan(4),atnami,atnamj
      character*8 tname,resnam,rn
      character*200 inpfile,analfile
      character* 132 line(maxrec)
      character*2 iatnm2
      common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
     -  mmatno(64),iatnm2(99)
      character*11 trajformatname
      common /trajectory/ nmmccheck,iftrajtyp(6),trajformatname(6)
      character*200 trajnam,trajnam2
      common /trajname/ trajnam,trajnam2,ltrajnam,ltrajnam2,ifirsttraj,
     -  ifirsttraj2,ilasttraj,ilasttraj2,incrementtraj,incrementtraj2
      character*4 namfcg
      character*4 tanames
      character*8 tnames
      common /tordat/ ntorn,tanames(4,28),tnames(28)
      common /connatdat/ ramax(99),ramax2(99),hlimfac,ianfg(99),
     -  namfcg(100),nrmw
      common /graphics/ npixx,npixy,maxpixx,maxpixy,idwmain,idwplot,
     -  wx,wy,wz,wxdr
      common /rotmat/ matrot0(4,4),matrot(4,4),nomat0
      common /columnlim/ incol(19),iidcol(19),iialtcol(19),iiinscol(19),
     -  iinamcol(2,19),iirescol(2,19),iiccol(2,19),iiresncol(2,19),
     -  iiseqncol(2,19),iisegcol(2,19),iiresidcol(2,19),iiqcol(2,19),
     -  iipotcol(2,19),iiocccol(2,19),iichemcol(2,19)
      character*5 crdext
      common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
     -  iommod,iommc,iommc4,iogro,iomol2,iomae,iocif,ioxxx,ioins,
     -  ionxyz,iosxyz,iosxyzrq,iograsp,iofull,lext(19),crdext(19)
      parameter (MAXHX=50)
      common /prokink/ icab(MAXHX),icaa(MAXHX),icb(MAXHX),ica(MAXHX),
     -  inb(MAXHX),ina(MAXHX),icapr,icpr,inpr,nra,nrb,icbpr,icgpr,icdpr,
     -  iprintpk
c     All arrays are of length maxframe
      parameter (MAXFRAMES=50000,MAXCOPY=600)
      parameter (MAXCOPY6=MAXCOPY-6)
      common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
     -  increment,increment2,bdxy(2,MAXFRAMES),wbxy(2,MAXFRAMES),
     -  fsxy(2,MAXFRAMES),wfsxy(2,MAXFRAMES),psxy(2,MAXFRAMES),
     -  turnpr(2,MAXFRAMES),scalardat(2,MAXFRAMES,MAXCOPY6),
     -  x0(MAXCOPY),y0(MAXCOPY),nxselres,ixselres(MAXCOPY)
      nnamcol=inamcol2-inamcol1+1
      if (nconfig .eq. 1) then
        call quiz(ansrun,nntyp,' ',' ',0,'Bond information source',23,0,
     -    5,6,00)
        if (nntyp .gt. 1 ) then
          call top_to_bond(nntyp,nneig,nhneig,ineig,iatnum,n,0,
     -      itemp1,itemp2,maxng,maxrec)
        else
          print *,'Bonds are defined from the coordinates read'
          call nnlist(nslt,islvw,naslv,n,iatnum,ifchrg,c,nneig,nneiga,
     -      nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,npneig,line,
     -      irescol1,irescol2,inamcol1,inamcol2,index,nconfig,innlist,
     -      molresflag,hblimfac,angmin,0,ibnd,indexo,isegno,ixres,
     -      maxrepconf,0,0,radtodeg,0,maxbox,maxng,maxrsd,maxrec)
        end if
        namleno=0
        call quiz(ansrun,iansrun,' ','syntax for',10,
     -    'torsion input generation',24,0,5,6,0)
        if (ansrun .eq. 't') then
          itform=2
          extnam='.mmt'
        else if (ansrun .eq. 'p') then
          itform=3
          extnam='.mmp'
        else if (ansrun .eq. 'r') then
          itform=4
          extnam='.mcs'
        else if (ansrun .eq. 'i') then
          itform=5
          extnam='.mcp'
        else if (ansrun .eq. 'g') then
          itform=6
          extnam='.mcg'
        end if
c       if (extnam(2:3) .eq. 'mc') call askyn(
c    -      'Do you want to include torsions moving hydrogens only',53,
c    -      0,-1,nohrot,0,0)
        if (extnam(2:3) .eq. 'mc')
     -    call getint('Maximum number of terminal hydrogens to rotate',
     -      46,0,1,3,maxhrot,133)
        analfile=inpfile
        if (namleno .ne. 0) close (40)
        namleno=namleni
        analfile(namleno-3:namleno)=extnam
        call openfile(40,0,'analysis',8,'new',analfile,namleno,
     -    notfnd,0,1,1,0,0)
      else
        write (40,2055) nconfig
      end if
      ic1=inamcol1
      ic2=inamcol2
      ir1=iirescol(1,inpcrdtyp)
      ir2=iirescol(2,inpcrdtyp)
      nrescol=ir2-ir1+1
      if (itform .ge. 4 .and. itform .le. 6) then
c       MMC protein torsion input
        if (nconfig .eq. 1) write (6,2001) 'MMC',analfile(1:namleno)
        ntorwr=0
        do ia=1,n
c         write (77,9781) ia,iatnum(ia),nneig(ia),nhneig(ia),mask(ia)
c9781     format(i5,' iatno=',i2,' nn=',i2,' nnh=',i2,' mask=',i2)
          if (iatnum(ia) .gt. 1 .and. nneig(ia) .gt. 1 .and.
     -      (nneig(ia)-nhneig(ia) .gt. 1 .or. nhneig(ia) .le. maxhrot)
     -      .and. mask(ia) .gt. 0) then
            resnam(1:nrescol)=line(index(ia))(ir1:ir2)
            call leftadjust4(line(index(ia))(ic1:ic2),atnami)
            do j=1,nneig(ia)
              ja=ineig(j,ia)
c             nhneig0=nhneig(ja)*nohrot
c             if (ia .lt. ja .and. nneig(ia)-nhneig(ia) .gt. 1 .and.
c    -            nneig(ja)-nhneig0 .gt. 1) then
              if (iatnum(ja) .gt. 1 .and. nneig(ja) .gt. 1 .and.
     -           ia .lt. ja .and. (nneig(ja)-nhneig(ja) .gt. 1 .or.
     -             nhneig(ja) .le. maxhrot) .and. mask(ja) .gt. 0) then
                call leftadjust4(line(index(ja))(ic1:ic2),atnamj)
c               Bond ia-ja found (on a heavy atom chain when at most maxhrot
c               hydrogens are moved)
                call checktorbond(resnam,ixres(ia),ixres(ja),atnami,
     -            atnamj,fixbond,ipep,ican,icac,issb)
c             write (77,9877) ia,ixres(ia),resnam,ja,ixres(ja),atnamj,
c    -           fixbond
c9877         format(i5,1x,i5,1x,a,' - ',i5,1x,i5,1x,a,' FIXBOND=',i2)
c               Loop and normal torsion stepsizes
                if (fixbond .eq. 0) then
                  dloop=30.0
                  dtor=30.0
                  if (itform .eq. 6 .or.
     -              (itform .eq. 4 .and. ipep+ican+icac+issb .eq. 0)
     -              .or. (itform .eq. 5 .and. ipep .eq. 0)) then
                    ntorwr=ntorwr+1
                    call readint(line(index(ia)),iresncol1,iresncol2,
     -                irna,2,1,irerr)
                    write (40,2051) ia,ja,dloop,dtor,ntorwr,
     -                irna,line(index(ia))(ir1:ir2),atnami,
     -                line(index(ja))(ir1:ir2),atnamj
                  end if
                end if
              end if
            end do
          end if
        end do
        if (nconfig .le. maxrepconf)
     -    print *,'Number of torsion angles written out=',ntorwr
      else
c       Macromodel torsion input
        if (nconfig .eq. 1)
     -    write (6,2001) 'Macromodel',analfile(1:namleno)
        ntorwr=0
        do ia=1,nslt
          do j=1,nneig(ia)
            ja=ineig(j,ia)
            if (ia .lt. ja .and. mask(ia)*mask(ja) .gt. 0) then
c             Bond ia-ja found
              do ii=1,nneig(ia)
                iaa=ineig(ii,ia)
                do jj=1,nneig(ja)
                  jaa=ineig(jj,ja)
                  if (jaa .ne. ia .and. iaa .ne. ja
     -                .and. iaa .ne. jaa) then
                    tors=dihangl(c,iaa,ia,ja,jaa,1,maxrec)*(180.0/pi)
                    tname='xxxxx'
                    if (inpcrdtyp .le. ioins) then
                      call leftadjust4(line(index(iaa))(ic1:ic2),lan(1))
                      call leftadjust4(line(index(ia))(ic1:ic2),lan(2))
                      call leftadjust4(line(index(ja))(ic1:ic2),lan(3))
                      call leftadjust4(line(index(jaa))(ic1:ic2),lan(4))
                      do i=1,ntorn
                        do k=1,4
                          if (lan(k) .ne. tanames(k,i)) go to 9011
                        end do
                        tname=tnames(i)
                        go to 9013
9011                    do k=1,4
                          if (lan(4-k+1) .ne. tanames(k,i)) go to 9012
                        end do
                        tname=tnames(i)
                        go to 9013
9012                    continue
                      end do
9013                  if (itform .eq. 3) then
c                       Protein, eliminate xxxxx, omega, and PRO phi,chi1,chi2
                        if (tname .eq. 'xxxxx' .or. tname .eq. 'omega')
     -                    then
                          tname='        '
                        else
                          resnam='        '
                          resnam(1:nrescol)=line(index(ia))(ir1:ir2)
                          call leftadjustn(resnam,rn,8)
                          if (rn .eq. 'PRO     ' .and.
     -                      (tname .eq. 'phi      ' .or.
     -                       tname .eq. 'chi1    '
     -                      .or. tname .eq. 'chi2    '))tname='        '
                        end if
                      end if
                      if (tname .ne. '        ') then
                        write (40,2053) iaa,ia,ja,jaa,tname,
     -                    line(index(iaa))(ic1:ic2),
     -                    line(index(iaa))(ir1:ir2),
     -                    line(index(ia))(ic1:ic2),
     -                    line(index(ia))(ir1:ir2),
     -                    line(index(ja))(ic1:ic2),
     -                    line(index(ja))(ir1:ir2),
     -                    line(index(jaa))(ic1:ic2),
     -                    line(index(jaa))(ir1:ir2),tors
                        ntorwr=ntorwr+1
                      end if
                    else
                      write (40,2054) iaa,ia,ja,jaa,tname
                      ntorwr=ntorwr+1
                    end if
                  end if
                end do
              end do
            end if
          end do
        end do
        if (nconfig .le. maxrepconf)
     -    print *,'Number of torsion angles written out=',ntorwr
        if (itform .eq. 3)
     -    print *,'Examine carefully if all the C and N terminal ',
     -    'torsions are present'
      end if
      return
1000  format(a)
2001  format(1x,a,' torsion input is written to ',a)
2051  format(2i5,2f10.4,i5,'    1       1.0     180.0 !',
     -  ' iresa=',i5,1x,a4,1x,a4,'-',a4,1x,a4)
2053  format(' ITOR',i8,3i7,'     0.0000   180.0000',23x,'!',a5,1x,
     -       3('(',a4,1x,a4,')-'),'(',a4,1x,a4,')',
     -  ' ta=',f7.1)
2054  format(' ITOR',i8,3i7,45x,a4)
2055  format(1x,79('-'),/,' Torsion input from configuration # ',i5)
      end
