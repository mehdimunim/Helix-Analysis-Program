      subroutine createrec(line,inpcrdtyp,ioutyp,cx,cy,cz,altcol,inscol,
     -  atnam,resnam,segnam,iseqno,iresnum,iresid,chemnam,potnam,
     -  frocc,q,nqdec,iqspace,mmcgm,nn,in,nconfig,ibnd,ihetat,blankline)
      character* 132 line,blankline
      character*1 altcol,inscol,resnam1
      character*4 atnam,segnam,chemnam
      character*6 potnam
      character*8 resnam
      character*5 crdext
      character*2 mmcgm
      common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
     -  iommod,iommc,iommc4,iogro,iomol2,iomae,iocif,ioxxx,ioins,
     -  ionxyz,iosxyz,iosxyzrq,iograsp,iofull,lext(19),crdext(19)
      dimension in(nn),ibnd(nn)
c      write (77,7711) ioutyp,potnam,atnam,resnam,iresnum
c7711  format(' CREATEREC ioutyp,potnam,atnam,resnam=',i4,3(2x,a5),
c    -   ' iresnum=',i5)
c      write (77,*) 'CREATEREC ioutyp,iogro=',ioutyp,iogro
c      write (80,7611) iresnum,iresid,q,segnam,altcol
c7611  format(' CREATEREC ires,resid=',2i5,' q=',f10.5,' segnam=',a,
c    -   ' ALTCOL=',a,'|')
c     print *,'CREATEREC inpcrdtyp,ioutyp=',inpcrdtyp,ioutyp
      call setcol(inpcrdtyp,ncol,idcol,ialtcol,iinscol,
     -  namcol1,namcol2,irescol1,irescol2,iccol1,iccol2,
     -  iresncol1,iresncol2,iseqncol1,iseqncol2,isegcol1,isegcol2,
     -  iresidcol1,iresidcol2,iqcol1,iqcol2,ipotcol1,ipotcol2,
     -  iocccol1,iocccol2,ichemcol1,ichemcol2,nrescol_i,nresncol,
     -  nsegcol_i,nnamcol_i,iofull)
      call setcol(ioutyp,ncol,idcol,ialtcol,iinscol,
     -  namcol1,namcol2,irescol1,irescol2,iccol1,iccol2,
     -  iresncol1,iresncol2,iseqncol1,iseqncol2,isegcol1,isegcol2,
     -  iresidcol1,iresidcol2,iqcol1,iqcol2,ipotcol1,ipotcol2,
     -  iocccol1,iocccol2,ichemcol1,ichemcol2,nrescol,nresncol,
     -  nsegcol,nnamcol,iofull)
      line=blankline
      if (ioutyp .eq. iocha) then
        icinc=0
        if (icinc .gt. 999999) icinc=1
        write (line(iseqncol1+icinc:iseqncol2+icinc),1000)
     -    mod(iseqno,100000)
        write (line(iresncol1+icinc:iresncol2+icinc),1000)
     -    mod(iresnum,100000)
        if (iseqno .eq. 100000) write (6,1100) 'seqno',99999,'CRD'
        if (iresnum .eq. 100000) write (6,1100) 'resnum',99999,'CRD'
        write (line(iccol1+icinc:iccol2+icinc),1011) cx,cy,cz
        call putreal(line(iqcol1+icinc:iqcol2+icinc),iqcol2-iqcol1+1,
     -    q,nqdec)
        write (line(iresidcol1+icinc:iresidcol2+icinc),1001)
     -    mod(iresid,10000)
        if (iresid .eq. 10000) write (6,1100) 'iresid',9999,'Charmm'
        call leftadjustline(line,iresidcol1+icinc,iresidcol2+icinc)
      else if (ioutyp .eq. iochaex) then
        write (line(iseqncol1:iseqncol2),1005) iseqno
        write (line(iresncol1:iresncol2),1005) iresnum
        write (line(iccol1:iccol2),1015) cx,cy,cz
        call putreal(line(iqcol1:iqcol2),iqcol2-iqcol1+1,q,2*nqdec)
        write (line(iresidcol1:iresidcol2),1004) iresid
        call leftadjustline(line,iresidcol1,iresidcol2)
      else if (ispdb(ioutyp) .gt. 0) then
        if (resnam .eq. '     ' .or. resnam .eq. 'UNK  ') then
          line(1:6)='HETATM'
        else
          line(1:4)='ATOM'
          if (ihetat .eq. 1) line(1:6)='HETATM'
        end if
        write (line(iseqncol1:iseqncol2),1000) mod(iseqno,100000)
        write (line(iresncol1:iresncol2),1001) mod(iresnum,10000)
        write (line(iccol1:iccol2),1012) cx,cy,cz
        line(ialtcol:ialtcol)=altcol
        line(iinscol:iinscol)=inscol
        if (nconfig .le. 10) then
          if (iresnum .eq. 10000) write (6,1100) 'resnum',9999,'PDB'
          if (iseqno .eq. 100000) write (6,1100) 'seqno',99999,'PDB'
        end if
        if (frocc .eq. 1.0) then
          line(55:60)='   1.0'
        else
          call putreal(line(55:60),6,frocc,0)
        end if
        if (iqspace .eq. 1) then
          line(iqcol1:iqcol1)=' '
          call putreal(line(iqcol1+1:iqcol2),iqcol2-iqcol1,q,nqdec)
        else
          call putreal(line(iqcol1:iqcol2),iqcol2-iqcol1+1,q,nqdec)
        end if
        if (ioutyp .eq. iocpdb) then
c         Right shift names that are less than four characters
          if (atnam(1:1) .ne. ' ' .and. atnam (4:4) .eq. ' ') then
            atnam(2:4)=atnam(1:3)
            atnam(1:1)=' '
          end if
          if (ichemcol2 .gt. ichemcol1)
     -      call writeitem(line,ichemcol1,ichemcol2,chemnam,0)
        end if
      else if (ioutyp .eq. iommod) then
        write (line(iresncol1:iresncol2),1000) iresnum
        write (line(iccol1:iccol2),1013) cx,cy,cz
        call writeitem(line,ipotcol1,ipotcol2,potnam,0)
        line(110:118)='  0.00000'
        line(iqcol1:iqcol2)=' '
        call putreal(line(iqcol1+1:iqcol2),iqcol2-iqcol1,q,nqdec)
        do i=1,nn
          write (line(12+(i-1)*8:12+(i-1)*8),1002) ibnd(i)
          write (line(5+(i-1)*8:10+(i-1)*8),1003) in(i)
        end do
        do i=nn+1,6
          line(4+(i-1)*8+1:4+i*8)='     0 0'
        end do
c       Create also the 1-letter residue code
        call changeprot(resnam,resnam1,2)
        line(isegcol1-1:isegcol2-1)=resnam1
      else if (ioutyp .eq. iommc) then
        write (line(iresncol1:iresncol2),1000) iresnum
        write (line(iccol1:iccol2),1011) cx,cy,cz
        call writeitem(line,ipotcol1,ipotcol2,potnam,0)
        call putreal(line(iqcol1:iqcol2),iqcol2-iqcol1+1,q,nqdec)
        line(51:52)=mmcgm
      else if (ioutyp .eq. iogro) then
        write (line(iseqncol1:iseqncol2),1000) iseqno
        write (line(iresncol1:iresncol2),1000) iresnum
        line(irescol1:irescol1)=' '
        line(namcol1:namcol1)=' '
        write (line(iccol1:iccol2),1012) cx/10.0,cy/10.0,cz/10.0
        write (line(iccol1+24:iccol2+24),1012) 0.0,0.0,0.0
      else if (ioutyp .eq. ioins) then
        write (line(iresncol1:iresncol2),1001) iresnum
        call writeitem(line,ichemcol1,ichemcol2,chemnam,0)
        if (inpcrdtyp .eq. ioutyp) then
          call writeitem(line,ipotcol1,ipotcol2,potnam,0)
        else
          call writeitem(line,ipotcol1,ipotcol2,chemnam,0)
        end if
        write (line(iccol1:iccol2),1014) cx,cy,cz
        call putreal(line(iqcol1:iqcol2),iqcol2-iqcol1+1,q,nqdec)
      end if
      call writeitem(line,irescol1,irescol2,resnam,nrescol_i)
      call writeitem(line,isegcol1,isegcol2,segnam,min0(4,nsegcol_i))
      call writeitem(line,namcol1,namcol2,atnam,min0(4,nnamcol_i))
      return
1000  format(i5)
1001  format(i4)
1002  format(i1)
1003  format(i6)
1004  format(i8)
1005  format(i10)
1011  format(3f10.5)
1012  format(3f8.3)
1013  format(3f12.5)
1014  format(3f15.9)
1015  format(3f20.10)
1100  format(' WARNING: leading ',a,' digits over ',i7,
     -  '  are dropped (',a,' format limit)')
      end
