      subroutine gcgwrite(inpcrdtyp,line,nhead,inpfile,leninpfile,ifile,
     -  seq,lenseq,skipnuc,maxrec)
      character*1 seq(lenseq),type
      character*4 pdbid
      character* 132 line(maxrec)
      character*200 inpfile
      character*40 keywords
      character*62 textinp
      character*5 crdext
      common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
     -  iommod,iommc,iommc4,iogro,iomol2,iomae,iocif,ioxxx,ioins,
     -  ionxyz,iosxyz,iosxyzrq,iograsp,iofull,lext(19),crdext(19)
      logical skipnuc
c     print *,'GCGWRITE lenseq=',lenseq,' ifile=',ifile
c     Determine if protein or nucleic acid
      type='N'
      do i=1,lenseq
        if (seq(i) .ne. 'A' .and. seq(i) .ne. 'C' .and.
     -      seq(i) .ne. 'G' .and. seq(i) .ne. 'T' .and.
     -      seq(i) .ne. 'U') then
c         Amino acid residue found
          type='P'
          go to 31
        end if
      end do
31    if (skipnuc .and. type .eq. 'N') return
      pdbid='    '
      if (inpcrdtyp .gt. 0) then
        if (ispdb(inpcrdtyp) .gt. 0) then
c         Translate header into keywords and pdb id
          nhline=0
          do i=1,nhead
            if (line(i)(1:6) .eq. 'HEADER') then
              nhline=nhline+1
              if (nhline .eq. 1) then
                pdbid=line(i)(63:66)
                keywords=line(i)(11:50)
              end if
            end if
          end do
        end if
      end if
      if (pdbid .eq. '    ') then
        call getname(keywords,len,'Keywords (max 40 chars)',23,40,'',0,
     -    0,0,0)
        call getname(pdbid,len,
     -    '4-character identifier (for LOCUS key)',38,4,'',0,0,0,0)
      end if
      write (ifile,1001) pdbid,lenseq
c     Calculate the checksum
      icheck=0
      icount=0
      do l=1,lenseq
        icount=icount+1
        lchar=ichar(seq(l))  !convert to ascii decimal form
        if (lchar.ge.97 .and. lchar.le.122)
     -    lchar=lchar-32     !change lower case to upper case
        icheck=icheck + icount*lchar
        if (icount.eq.57) icount=0
      end do
      lsum=icheck/10000                  !integer division
      icheck=icheck - lsum*10000        !remainder which is the checksum
      ncline=0
      if (inpcrdtyp .gt. 0) then
        if (ispdb(inpcrdtyp) .gt. 0) then
c         Get compound name as definition
          do i=1,nhead
            if (line(i)(1:6) .eq. 'COMPND') then
              ncline=ncline+1
              if (ncline .eq. 1) then
                write (ifile,1003) line(i)(11:72)
              else
                write (ifile,1004) line(i)(11:72)
              end if
            end if
          end do
        end if
      end if
      if (ncline .eq. 0) then
        call getname(textinp,len,'Compound name (for DEFINITION key)',
     -    34,62,'',0,0,0,3)
        write (ifile,1003) textinp
      end if
      write (ifile,1002) inpfile(1:leninpfile)
      write (ifile,1006) keywords
      write (ifile,1005) pdbid,lenseq,type,icheck
      nline=(lenseq-1)/50+1
      if=1
      linelim=50
      do iline=1,nline
        linc=(iline-1)*50
        if (iline .eq. nline) linelim=lenseq-linc
        nfirst=(iline-1)*50+1
        write (ifile,1000) nfirst,(seq(linc+l),l=if,linelim)
      end do
      return
1000  format(i8,1x,5(10a1,1x))
1001  format('LOCUS       ',a4,'------',i7,' AA    PROT',
     -  12x,'PRE-ENTRY 00/00/0')
1002  format('ORIGIN      ',a)
1003  format('DEFINITION  ',a62)
1004  format(12x,a62)
1005  format(1x,a4,' Length:',i5,' April 2, 1992      Type: ',a1,
     -  '  Check:',i6,' ..')
1006  format('KEYWORDS    ',a62)
      end
