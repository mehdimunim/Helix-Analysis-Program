      subroutine convertseq(resnames,resnames1,segnames,inpfile,outfile,
     -  maxrsd)
c    -  filename,namlen,notfound,idatapath,ib,nosys,noecho,ioverall)
      character*8 resnames(maxrsd)
      character*1 resnames1(maxrsd)
      character*200 inpfile,outfile
      character*4 segnames(maxrsd) 
      character*1 ans,ch
      character*4 s4
      character*132 l132(1)
      character*200 line
      character*1 abc,digits,hexdigits
      common /charactersets/ ihex(25),abc(62),digits(14),hexdigits(25)
      linpfile=0
      loutfile=0
      call openfile(30,0,'input sequence',14,'old',inpfile,linpfile,
     -  notfndi,0,1,1,0,0)
      call quiz(ans,inpstyp,' ',' ',0,
     -  'Input format for sequence list',30,0,5,6,00)
      call openfile(31,0,'output sequence',15,'new',outfile,loutfile,
     -  notfndo,0,1,1,0,1)
      call quiz(ans,iostyp,' ',' ',0,
     -  'Output format for sequence list',31,0,5,6,00)
      if (notfndi+notfndo .gt. 0) stop
      nres=0
      nresp=nres
      nseg=0
      if (inpstyp .eq. 1) then
c       Charmm sequence input
        print *,'It is assumed that no other Charmm input items are ',
     -    'in the file'
        do while (.true.)
          call blankout(line,1,80)
          read (30,1000,end=100) line
          call lastchar(line,lc,80)
          if (line(1:1) .eq. '*' .and. lc .eq. 1) then
c           Header read
            read (30,*,end=888,err=888) numres
            nr=0
            do while (nr .lt. numres)
              call blankout(line,1,80)
              read (30,1000,end=100) line
              call lastchar(line,lc,80)
              if (mod(lc,5) .ne. 0) then
                print *,'List of residues seem to be corrupted:' 
                print *,line(1:lc)
                stop
              end if
              nc=lc/5
              do i=1,nc
                nres=nres+1
                resnames(nres)(1:3)=line((i-1)*5+3:i*5)
                call changeprot(resnames(nres)(1:3),resnames1(nres),2)
              end do
              nr=nr+nc
            end do
            call blankout(line,1,80)
            read (30,1000,end=100) line
            nseg=nseg+1
            s4='    '
            if (line(1:8) .ne. 'GENERATE') then
c             print *,'GENERATE * SETUP record is missing'
              s4(1:1)=abc(nseg)
            else
              ic=10
              call nextblank(line,ic,80)
              len=min0(4,ic-10)
              s4(1:len)=line(10:9+len)
            end if
            do i=1,numres
              segnames(nresp+i)=s4
            end do
            nresp=nres
          end if
        end do
      else if (inpstyp .eq. 2) then
c       PDB SSEQRES input
        ch=' '
        print *,'SEQRES input start'
        do while (.true.)
          call blankout(line,1,80)
          read (30,1000,end=100) line
          if (line(1:6) .eq. 'SEQRES') then
            call lastchar(line,lc,80)
            if (line(12:12) .ne. ch) nseg=nseg+1
            ch=line(12:12)
c           print *,'L12=',line(12:12),' CH=',ch,' NS=',nseg,' NR=',nres
            ic0=20
            do while (ic0 .lt. lc)
              nres=nres+1
              segnames(nres)(1:1)=ch
              segnames(nres)(2:4)='   '
              resnames(nres)(1:3)=line(ic0:ic0+2)
              call changeprot(resnames(nres)(1:3),resnames1(nres),2)
              ic0=ic0+4
            end do
          end if
        end do
      else if (inpstyp .eq. 3 .or. inpstyp .eq. 4) then
c       >Title + <1>-char list (FASTA) or PIR
        nseg=1
        call blankout(line,1,80)
        do while (line(1:1) .ne. '>')
          call blankout(line,1,80)
          read (30,1000,end=100) line
        end do
        do while (.true.)
          if (line(1:6) .ne. 'SEQNAM') then
            call lastchar(line,lc,80)
            do ic=1,lc
              nres=nres+1
              resnames1(nres)=line(ic:ic)
              call changeprot(resnames(nres)(1:3),resnames1(nres),1)
              segnames(nres)='A   '
            end do
          end if
        end do
      else if (inpstyp .eq. 5) then
c       CIF
        do while (.true.)
          call blankout(line,1,200)
          read (30,1000,end=100) line
          call lastchar(line,lc,200)
          ic=1
          call findnextchar('p',line,ic,200)
          if (lc-11 .gt. ic) then
            if (line(ic:ic+10) .eq. 'polypeptide') then
              ch=abc(nseg+1)
              if (line(lc:lc) .eq. '?') then
c               Sequence in the same line
                do k=1,3
                  call nextblank(line,ic,200)
                  ic=ic+1
                end do
                call nextchar(line,ic,200)
                do while (line(ic:ic) .ne. ' ')
                  nres=nres+1
                  resnames1(nres)=line(ic:ic)
                  segnames(nres)=ch//'   '
                  call changeprot(resnames(nres)(1:3),resnames1(nres),1)
                  ic=ic+1
                end do
              else
                nl=0
                do while (nl .eq. 0 .or. line(1:1) .ne. ';')
                  call blankout(line,1,200)
                  read (30,1000,end=100) line
                  if (nl .eq. 0 .or. line(1:1) .ne. ';') then
                    call lastchar(line,lc,200)
                    incr=0
                    if (nl .eq. 0) incr=1
                    do ic=1+incr,lc
                      resnames1(nres-incr+ic)=line(ic:ic)
                      segnames(nres-incr+ic)=ch//'   '
                      call changeprot(resnames(nres-incr+ic)(1:3),
     -                  resnames1(nres-incr+ic),1)
                    end do
                    nres=nres+lc-incr
                    nl=nl+1
                    line(1:1)=' '
                  end if
                end do
              end if
              nseg=nseg+1
            end if
          end if
        end do
      else if (inpstyp .eq. 6) then
c       GCG
        s4='    '
        do while (.true.)
          call blankout(line,1,200)
          read (30,1000,end=100) line
          call lastchar(line,lc,200)
          if (line(1:5) .eq. 'LOCUS') then
            ic=6
            call nextchar(line,ic,200)
            call nextblank(line,ic,200)
            call nextchar(line,ic,200)
            ic1=ic
            call nextblank(line,ic,200)
            read (line(ic1:ic-1),*,err=888) numres
            do i=1,4
              read(30,1000,end=100) 
            end do
            nr=0
            do while (nr .lt. numres)
              ir1=nr+1
              ir2=min0(nr+50,numres)
              read (30,1000,end=100) line
              read(line,2001) (resnames1(nresp+i),i=ir1,ir2)
              nr=nr+50
            end do
            nseg=nseg+1
            s4(1:1)=abc(nseg)
            do i=nresp+1,nresp+numres
              segnames(i)=s4
              call changeprot(resnames(i)(1:3),resnames1(i),1)
            end do
            nresp=nresp+numres
            nres=nresp
c           print *,'NSEG=',nseg,' NRES=',nres,' S4=',s4
          end if
        end do
      end if
100   print *,'Read ',nres,' residues in ',nseg,' segments'
      nresp=0
      segnames(nres+1)='    '
      do iseg=1,nseg
        s4=segnames(nresp+1)
        nrespc=nresp
        do while (segnames(nrespc+1) .eq. s4 .and. nrespc .lt. nres)
          nrespc=nrespc+1
        end do 
c       print *,'ISEG=',iseg,' INPSTYP=',inpstyp,' IOSTYP=',iostyp
        if (iostyp .eq. 1) then
c         Charmm sequence format
          write (31,2024) s4,nrespc-nresp,
     -      (resnames(i)(1:3),i=nresp+1,nrespc)
          write (31,2021) s4
        else if (iostyp .eq. 2) then
c         PDB SSEQRES iout
          ir1=nresp+1
          nl=0
          do while (ir1 .le. nrespc)
            nl=nl+1
            ir2=min0(nrespc,ir1+12)
            write (31,2100) nl,s4(1:1),nrespc-nresp,
     -        (resnames(ir)(1:3),ir=ir1,ir2)
            ir1=ir2+1
          end do
        else if (iostyp .eq. 3) then
c         >Title + <1>-char list (FASTA)
          if (iseg .eq. 1) write (31,1000) '>SIMULAID conversion'
          write (31,2000) (resnames1(i),i=nresp+1,nrespc)
        else if (iostyp .eq. 4) then
c         PIR
          if (iseg .eq. 1) write (31,1000) '>SIMULAID conversion'
          write (31,2000) (resnames1(i),i=nresp+1,nrespc)
        else if (iostyp .eq. 5) then
c         GCG
          if (inpstyp .eq. 2) write (6,2029)
          call gcgwrite(0,l132,0,outfile,loutfile,31,
     -      resnames1(nresp+1),nrespc-nresp,.false.,1)
        end if
        nresp=nrespc
      end do
      close(30)
      close(31)
      return
888   write (6,*) 'Invalid number of residues ic=',ic,' in ',line(1:40)
      stop
1000  format(a)
2000  format(75a1)
2001  format(9x,10a1,1x,10a1,1x,10a1,1x,10a1,1x,10a1,1x)
2021  format('GENERATE ',a4,' SETUP',/)
2024  format('* Sequence list from SIMULAID Segid=',a4,/,
     -  '* ',/,i5,/,(12a5))
2029  format(' NOTE: If the input is a full PDB file than the ',
     -  '<E>xtract sequence',/,' option can also obtain the additional',
     -  ' information asked below')
2100  format('SEQRES',i4,1x,a1,i5,1x,13(1x,a3))
      end
