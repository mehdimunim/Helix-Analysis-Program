      subroutine writeseq(outfile,namleno,numres,resnames,resnames1,
     -  numslv,line,index,inpcrdtyp,title,nslt,n,ionly,maxrsd,
     -  maxrec)
c*****Generates a residue list in a variety of formats
      dimension index(maxrec)
      character*200 outfile,seqfile
      character* 132 line(maxrec)
      character*80 title
      character*1 ans,resnames1(maxrsd)
      character*4 segid
      character*8 resnames(maxrsd)
      character*1 aanames1
      character*2 mmodtoamb
      character*3 aanames3
      character*80 seqid
      common /atnamcon/ mmodtoamb(100),aanames1(58),aanames3(58),
     -  naanames,nnanames,nnammnames,nnames,ixwatnam
      common /columnlim/ incol(19),iidcol(19),iialtcol(19),iiinscol(19),
     -  iinamcol(2,19),iirescol(2,19),iiccol(2,19),iiresncol(2,19),
     -  iiseqncol(2,19),iisegcol(2,19),iiresidcol(2,19),iiqcol(2,19),
     -  iipotcol(2,19),iiocccol(2,19),iichemcol(2,19)
      character*1 abc,digits,hexdigits
      common /charactersets/ ihex(25),abc(62),digits(14),hexdigits(25)
      character*5 crdext
      common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
     -  iommod,iommc,iommc4,iogro,iomol2,iomae,iocif,ioxxx,ioins,
     -  ionxyz,iosxyz,iosxyzrq,iograsp,iofull,lext(19),crdext(19)
      character*6 typnam(5)
      character*30 question
      data typnam /'Charmm','PDB   ','1-char','PIR   ','GCG   '/
c     print *,'WRITESEQ numres,numslv=',numres,numslv
      if (inpcrdtyp .gt. ioins) then
        print *,'This input format does not have sequence information'
        stop
      end if
      if (ionly .eq. 0) then
c       Run incidental to sequence generation - use default format w/o/ asking
        call askyn('Do you want a sequence list',27,1,-1,iseqls,0,0)
        if (iseqls .eq. 0) return
        iostyp=inpcrdtyp
        if (iostyp .gt. 2) iostyp=2
        write (6,2002) typnam(iostyp)
      else
c       Ask for format
        call quiz(ans,iostyp,' ',' ',0,
     -    'Output format for sequence list',31,0,5,6,00)
      end if
      call changeext(outfile,seqfile,namleno,namlens,'seq',3,0,0)
      call openfile(30,0,' ',1,'new',seqfile,namlens,notfnd,0,1,1,0,0)
      isegidc1=iisegcol(1,inpcrdtyp)
      isegidc2=iisegcol(2,inpcrdtyp)
      iresnc1=iiresncol(1,inpcrdtyp)
      iresnc2=iiresncol(2,inpcrdtyp)
      iresc1=iirescol(1,inpcrdtyp)
      iresc2=iirescol(2,inpcrdtyp)
      nrescol=iresc2-iresc1+1
c     Write sequence input (by segments)
      if (nslt .eq. n) then
        index(n+1)=index(n)+1
        line(index(n+1))(isegidc1:isegidc1)='*'
        line(index(n+1))(iresnc1:iresnc1)='*'
      end if
      lenseq=1
      lenseqo=1
      ia=1
      nseg=0
      nfail=0
      do while (ia .le. nslt)
        nseg=nseg+1
        segid=line(index(ia))(isegidc1:isegidc2)
        do while (segid .eq.
     -    line(index(ia))(isegidc1:isegidc2) .and. ia .le. nslt)
            ia=ia+1
            if (line(index(ia))(iresnc1:iresnc2) .ne.
     -          line(index(ia-1))(iresnc1:iresnc2)) then
              resnames(lenseq)(1:nrescol)=
     -          line(index(ia-1))(iresc1:iresc2)
              lenseq=lenseq+1
            end if
        end do
        if (iostyp .eq. 1) then
          write (30,2022)
          ifc=1
          call nextchar(title,ifc,80)
          if (ifc .lt. 78) write (30,2020) '* ', title(1:78)
          write (30,2024) segid,lenseq-lenseqo,
     -      (resnames(i)(1:nrescol),i=lenseqo,lenseq-1)
          write (30,2021) segid
          if (numslv .gt. 0) then
            write (30,2023)
     -        line(index(nslt+1))(isegidc1:isegidc2),numslv
            if (line(index(nslt+1))(isegidc1:isegidc2) .eq.
     -          line(index(nslt))(isegidc1:isegidc2))
     -      write (6,2026) line(index(nslt))(isegidc1:isegidc2)
          end if
        else if (iostyp .eq. 2) then
          ir1=lenseqo
          nl=0
          do while (ir1 .le. lenseq-1)
            nl=nl+1
            ir2=min0(lenseq-1,ir1+12)
            write (30,2100) nl,segid(1:1),lenseq-lenseqo,
     -        (resnames(ir)(1:3),ir=ir1,ir2)
            ir1=ir2+1
          end do
        else if (iostyp .ge. 3) then
c         Convert to 1-character ID-s
          do ir=lenseqo,lenseq-1
            call changeprot(resnames(ir),resnames1(ir),2)
            if (resnames1(ir) .eq. '*') nfail=nfail+1
          end do
          if (iostyp .eq. 3) then
            if (nseg .eq. 1) write (30,2028) title
            write (30,2001) (resnames1(i),i=lenseqo,lenseq-1)
          else if (iostyp .eq. 4) then
c           PIR
            question='Sequence id of segment   '
            write (question(23:25),1002) nseg
            call getname(seqid,namlen,question,25,80,abc(nseg),1,0,0,0)
            call lastchar(title,lentit,80)
            lentit=min0(80,lentit+4+namlen)
            write (30,2027) seqid(1:namlen),title(1:lentit),
     -        seqid(1:namlen)
            resnames1(lenseq)='*'
            write (30,2001) (resnames1(i),i=lenseqo,lenseq)
          else
            call gcgwrite(inpcrdtyp,line,index(1),seqfile,namlens,30,
     -        resnames1(lenseqo),lenseq-lenseqo,.false.,maxrec)
          end if
        end if
        lenseqo=lenseq
      end do
      close (30)
      write (6,2025) typnam(iostyp),numres-numslv,seqfile(1:namlens)
      if (numslv .gt. 0) print *,'Number of solvent residues=',numslv
      if (nfail .gt. 0) print *,'WARNING',nfail,' residues had no ',
     -  '1-character version'
      return
1002  format(i3)
2001  format(80a1)
2002  format(' Sequence is saved in ',a,' format',/,' Use the ',
     -  '<E>xtract sequence option to get it in other formats')
2020  format(a,a)
2021  format('GENERATE ',a4,' SETUP',/)
2022  format('READ SEQUENCE CARD')
2023  format('!Solvent segment',/,'READ SEQU',a4,1x,i6,/)
2024  format('* Sequence list from SIMULAID Segid=',a4,/,
     -  '* ',/,i5,/,(12a5))
2025  format(1x,a,' sequence (input) list of ',i6,' residues were',
     -  ' written to file ',/,5x,a)
2026  format(' WARNING: Solvent segmentid (',a4,') is the same as',
     -  ' the last solute segment id')
2027  format('>P1;',a,1x,a,/,'sequence:',a,':::::::0.0: 0.0')
2028  format('>',a,a)
2100  format('SEQRES',i4,1x,a1,i5,1x,13(1x,a3))
      end
