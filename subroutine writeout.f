      subroutine writeout(iout,inpcrdtyp,iotyp,line,index,isegno,n,
     -  marker,iwhead,imodelnum,ntitlin,ntitlinw,title,
     -  blankline,nosort,ianaltyp,ifromtraj,etot,ietot,nclstmem,noend,
     -  keeprem,iwriteatsym,iatnum,maxrec)
c     Puts out the generated configuration, with proper headers, separators
      dimension index(n),isegno(n),iatnum(n)
      character*6 marker(16)
      character*80 title,linew
      character* 132 line(maxrec),blankline,ansline,pline
      character*5 crdext
      common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
     -  iommod,iommc,iommc4,iogro,iomol2,iomae,iocif,ioxxx,ioins,
     -  ionxyz,iosxyz,iosxyzrq,iograsp,iofull,lext(19),crdext(19)
      common /logging/ logfile,ipredict
      character*11 formatname
      common /formats/ iqdconv(20),formatname(19)
      common /columnlim/ incol(19),iidcol(19),iialtcol(19),iiinscol(19),
     -  iinamcol(2,19),iirescol(2,19),iiccol(2,19),iiresncol(2,19),
     -  iiseqncol(2,19),iisegcol(2,19),iiresidcol(2,19),iiqcol(2,19),
     -  iipotcol(2,19),iiocccol(2,19),iichemcol(2,19)
      real*8 xtlabc,xtlabc0
      common /boxdat/ xtlabc(6),xtlabc0(6),box(3),box0(3),edge_gen(3,3),
     -  cell0(3,27),cell(3,27),cellalt(3,27),
     -  ncell,ioppbc,noboxinfoar,noboxinfow,noboxrep,
     -  istuner,iboxtypfound,ixcrd(3),ixang(3),ixyzhex(3),
     -  ixyzhextraj(3),isizewarn
      common /analparm/ nsltref_f,nsltref_l,rcut_cv,icvtyp
c     write (6,*) 'WRITEOUT inpcrdtyp,iotyp,n,iout=',
c    -  inpcrdtyp,iotyp,n,iout
      nline_add=0
      pline=blankline
      pline(1:6)=marker(iotyp)
      modelfound=0
      if (iwhead .eq. 1) then
        if (marker(iotyp) .ne. '      ') then
c         Write headers
          if (ispdb(iotyp) .gt. 0 .or. ischarmm(iotyp) .eq. 1) then
c           Charmm or PDB outputs
            if (inpcrdtyp .eq. iommod) then
              icol=60
              call nextblank(line(1),icol,132)
              pline(8:icol+1)=line(1)(7:icol)
              call writeline(iout,pline,1,icol+1,0)
              pline(icol+2:132)=blankline(icol+2:132)
              icol=icol-50
              call nextchar(line(1)(51:130),icol,132)
              call lastchar(line(1),ifc,130)
              if (icol .lt. 80) then
                pline(8:ifc-(icol+50)+8)=line(1)(icol+50:ifc)
                call writeline(iout,pline,1,ifc-(icol+50)+8,0)
              end if
            else
              if (ispdb(inpcrdtyp) .gt. 0) then
                mncl=7
                mxcl=79
              else
                mncl=2
                mxcl=74
              end if
              ilw=ntitlinw
              if (ntitlin .eq. 0) then
                call blankout(linew,1,80)
                if (ischarmm(iotyp) .eq. 1) then
                  write (linew,1000) marker(iotyp)
                else
                  write (linew,1000) 'TITLE'
                end if
                call lastchar(linew,lc0,80)
                call lastchar(title,lct,80)
                ltitle=min0(lct,80-lc0-1)
                lc=lc0+ltitle+1
                linew(1:lc)=linew(1:lc0)//' '//title(1:ltitle)
                write (iout,1000) linew(1:lc)
                nline_add=1
              end if
              do i=1,ntitlin
                if (line(i)(1:5) .eq. 'MODEL') modelfound=1
                if (ischarmm(iotyp) .eq. 1 .and. i .le. ntitlin) then
c                 Make sure not to write more than 30 lines, and no blanks
                  if (ianaltyp .eq. 13 .or. ianaltyp .eq. 20 .or.
     -                ianaltyp .eq. 22) then
                    nline_add=1
                  else if (ianaltyp .eq. 14) then
                    nline_add=3
                  else
                    nline_add=0
                  end if
                  if (i .eq. 1 .and. (ietot .eq. 1 .or. nclstmem .gt. 0
     -                .or. ifromtraj .gt. 0)) then
                    nline_add=nline_add+1
                    call blankout(linew,1,80)
                    write (linew,1000) marker(iotyp)
                    call lastchar(linew,lc,80)
                    if (ietot .eq. 1) then
                      write (linew,2034) '*',etot
                      call lastchar(linew,lc,80)
                    end if
                    if (nclstmem .gt. 0) then
                      write (linew(lc+1:lc+16),2036) nclstmem
                      lc=lc+16
                    end if
                    if (ifromtraj .gt. 0) then
                      write (linew(lc+1:lc+16),2037) ifromtraj
                      lc=lc+16
                    end if
                    write (iout,1000) linew(1:lc)
                  end if
                  if (ilw .lt. 30-nline_add) then
c                   Skip energy record since that was only in the input template
                    ietotread=0
                    if (ifromtraj .gt. 0) call checkforetot(1,line(i),
     -                0,etotread,ietotread,0)
                    if (ietotread .eq. 1) nline_add=nline_add-1
                    if (ietotread .eq. 0) then
                      icol=mncl
                      call nextchar(line(i),icol,132)
                      call lastchar(line(i),ifc,mxcl)
                      if (icol .le. mxcl) then
                        ilw=ilw+1
                        pline(8:ifc-icol+8)=line(i)(icol:ifc)
                        call writeline(iout,pline,1,ifc-icol+8,0)
                      end if
                    end if
                  end if
                else
                  if (i .eq. 1 .and. (ietot .eq. 1 .or. nclstmem .gt. 0
     -                .or. ifromtraj .gt. 0)) then
                    nline_add=nline_add+1
                    call blankout(linew,1,80)
                    write (linew,1000) marker(iotyp)
                    call lastchar(linew,lc,80)
                    if (ietot .eq. 1) then
                      write (linew(lc+1:lc+19),2034) etot
                      lc=lc+19
                    end if
                    if (nclstmem .gt. 0) then
                      write (linew(lc+1:lc+16),2036) nclstmem
                      lc=lc+16
                    end if
                    if (ifromtraj .gt. 0) then
                      write (linew(lc+1:lc+16),2037) ifromtraj
                      lc=lc+16
                    end if
                    write (iout,1000) linew(1:lc)
                  end if
                  call lastchar(line(i),ifc,mxcl)
                  ietotread=0
                  if (ifromtraj .gt. 0)
     -              call checkforetot(6,line(i),0,etotread,ietotread,0)
                  if (ietotread .eq. 0) then
                    if (mncl .eq. 7) then
c                     PDB to PDB - keep full non-atom lines
                      call writeline(iout,line(i),1,ifc,0)
                    else
                      pline(8:ifc-mncl+8)=line(i)(mncl:ifc)
                      call writeline(iout,pline,1,ifc-mncl+8,0)
                    end if
                  end if
                end if
                pline(8:132)=blankline(8:132)
              end do
c             Add data column annotation, if required
              pline=blankline
              if (ianaltyp .eq. 13) then
                write (iout,2031) marker(iotyp),'hydropathy label'
              else if (ianaltyp .eq. 14) then
                if (icvtyp .eq. 1) then
                  write (iout,2031) marker(iotyp),'Circular variance'
                else
                  write (iout,2031) marker(iotyp),
     -              'weighted circular variance'
                end if
                write (pline,2032) marker(iotyp),rcut_cv
                call writeline(iout,pline,1,0,0)
                write (pline,2033) marker(iotyp),nsltref_f,nsltref_l
                call writeline(iout,pline,1,0,0)
              else if (ianaltyp .eq. 20) then
                write (iout,2031) marker(iotyp),
     -            'Delphi potential label'
              else if (ianaltyp .eq. 22) then
                write (iout,2031) marker(iotyp),
     -            'residue Root Mean Square Fluctuation (RMSF)'
              else if (ianaltyp .eq. 37) then
                write (iout,2031) marker(iotyp),
     -            'residue distance difference average'
              else if (ianaltyp .eq. 38) then
                write (iout,2031) marker(iotyp),
     -          'residue Root Mean Square Fluctuation (RMSF) difference'
              else if (ianaltyp .eq. 91) then
                write (iout,2035) marker(iotyp),'representative-based'
              else if (ianaltyp .eq. 92) then
                write (iout,2035) marker(iotyp),
     -             'closest approach-based'
              else if (ianaltyp .eq. 93) then
                write (iout,2035) marker(iotyp),'mutually proximal'
              else if (ianaltyp .eq. 94) then
                write (iout,2035) marker(iotyp),'sum of all'
              end if
            end if
            if (noboxinfoar .eq. 0) then
              write (iout,2030) marker(iotyp),box
            end if
          else if (iotyp .eq. ioins) then
            write (iout,2062)
            if (inpcrdtyp .le. ioins) then
              write (iout,1001) title
            else
              write (iout,2029) ' ',formatname(inpcrdtyp)
            end if
            write (iout,2063)
          end if
        end if
        if (inpcrdtyp .ne. iotyp .and. marker(iotyp) .ne. '      ')
     -    write (iout,2029) marker(iotyp),formatname(inpcrdtyp)
c       Write number of atoms
        if (iotyp .eq. iocha) then
          if (n .lt. 100000) then
            write (iout,1291) n
          else
            write (6,3000)
            call askstop(1)
            write (6,3001)
            write (iout,1292) n
          end if
        else if (iotyp .eq. iochaex) then
           write (iout,1293) n
        else if (iotyp .eq. iommod) then
          ansline=blankline
          write (ansline(1:6),1008) n
          if (title(1:4) .ne. '@#$%') then
            ansline(8:87)=title
          else
            ansline(8:28)='Converted by Simulaid'
            if (ntitlin . gt. 0)
     -        ansline(30:109)=title
          end if
          call lastchar(ansline,ifc,incol(inpcrdtyp))
          call writeline(iout,ansline,1,ifc,0)
        else if (iotyp .eq. iogro) then
          write (iout,1001) title
          write (iout,1007) n
        else if (iotyp .eq. ionxyz .or. iotyp .eq. iosxyz .or.
     -           iotyp .eq. iosxyzrq) then
          write (iout,1007) n
        end if
      end if
      if (ischarmm(iotyp) .eq. 1) then
c       Left-adjust residue ID
        iresidcol1=iiresidcol(1,iotyp)
        iresidcol2=iiresidcol(2,iotyp)
        do ia=1,n
          call leftadjustline(line(index(ia)),iresidcol1,iresidcol2)
        end do
      end if
      if (imodelnum .gt. 0) write (iout,1002) 'MODEL',imodelnum
      ncol=incol(iotyp)
      call lastchar(line(index(1)),ifc,ncol)
      if (n .gt. 0) then
        if (iwriteatsym .gt. 0)
     -    call addatsym(line(index(1)),iatnum(1),ifc)
        call writeline(iout,line(index(1)),1,ifc,0)
      end if
      if (nosort .eq. 1 .and.
     -    ispdb(inpcrdtyp) .gt. 0 .and. ispdb(iotyp) .gt. 0) then
c       REMARK records can only be placed safely if index is unchanged
        nrem=0
        do i=index(1),index(n)
          if (line(i)(1:6) .eq. 'REMARK') nrem=nrem+1
        end do
        if (nrem .gt. 0 .and. keeprem .eq. -1) then
          call askyn(
     -      'Do you want to keep REMARK records among the ATOM records',
     -      57,1,1,keeprem,000,0)
        end if
      end if
c      if (keeprem .gt. 0) write (77,8811)
c     -  (i,index(i),line(index(i))(1:80),i=1,index(n))
c8811  format(i5,' index-',i5,' line=',a)
c     ixprev=max0(0,index(1)-3)
      ixprev=index(1)
      do ia=2,n
        if (isegno(ia) .ne. isegno(ia-1)) then
c         Segment end
          if (ispdb(iotyp) .gt. 0) write (iout,1000) 'TER'
          if (iotyp .eq. ioins) write (iout,1000) 'end'
        end if
        if (keeprem .eq. 1) then
          if (index(ia)-ixprev .gt. 1) then
            do i=ixprev+1,index(ia)
              if (line(i)(1:6) .eq. 'REMARK') then
                call lastchar(line(i),lc,80)
                call writeline(iout,line(i),1,lc,0)
              end if
            end do
          end if
          ixprev=index(ia)
        end if
        call lastchar(line(index(ia)),ifc,ncol)
        if (iwriteatsym .gt. 0)
     -    call addatsym(line(index(ia)),iatnum(ia),ifc)
        call writeline(iout,line(index(ia)),1,ifc,0)
      end do
      if (imodelnum+modelfound .gt. 0) then
        write (iout,1000) 'ENDMDL'
      else if (ispdb(iotyp) .gt. 0 .and. noend .eq. 0) then
        write (iout,1000) 'END'
      end if
      if (iotyp .eq. ioins) write (iout,2064)
      return
1000  format(a)
1001  format(a80)
1002  format(a,i6)
1007  format(i5)
1008  format(i6)
1291  format('*',/,i5)
1292  format('*',/,i6)
1293  format('*',/,i10,2x,'EXT')
2029  format(a,' Simulaid generated this file from ',a,' input')
2030  format(a,' Cell dimensions: ',3f10.5)
2031  format(a,' Data column shows ',a)
2032  format(a,' Cutoff used for the CV calculation=',f5.1,' A')
2033  format(a,' Solute atom range used for the CV calculation [',i4,
     -  ',',i5,']')
2034  format(' E=',e16.8)
2035  format(a,' Data column shows ',a,' contact counts')
2036  format(' Nclstmem=',i6)
2037  format(' Frame #',i8)
2062  format('!BIOSYM archive 3',/,'PBC=OFF')
2063  format('!DATE')
2064  format('end',/,'end')
3000  format(' WARNING number of atoms exceeds the format limit',/,
     -  ' You may want to restart the run and ask for the extended ',
     -  ' Charmm format')
3001  format(' WARNING number of atoms will be printed with i6 format')
      end
