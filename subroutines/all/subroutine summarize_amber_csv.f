      subroutine summarize_amber_csv(resnames,minmaxlab1,mxrsd)
      character*8 resnames(mxrsd)
      character*1 minmaxlab1(mxrsd)
      character*1000 line
      parameter (MAXREC=200000,MAXRSD=70000,MAXPHI=400)
      parameter (MAXCOL=36)
      character*30 title(MAXCOL),label(4),syslabel
      character*30 source(3)
      character*1 ans,types1(4)
      character*26 err_typ(2)
      character*80 atitle
      character*200 inpstatfile,outtabfile
      parameter (IFILL7=MAXPHI*MAXPHI*MAXPHI-(2*MAXCOL+10)*MAXREC
     -  -6*MAXRSD)
      common /nnwork/ data(MAXCOL,MAXREC),err(MAXCOL,MAXREC),
     -  ixr1(MAXREC),ixr2(MAXREC),col1(MAXREC),col2(MAXREC),
     -  ixfres(MAXREC),ixlres(MAXREC),tabrowsum(MAXREC),ixcol(MAXREC),
     -  ixrow(MAXREC),iresix(MAXREC),colsum(MAXRSD),colsumm(MAXRSD),
     -  colmin(MAXRSD),colmax(MAXRSD),ixcolmin(MAXRSD),ixcolmax(MAXRSD),
     -  fill(IFILL7)
      dimension rmin(MAXCOL),rmax(MAXCOL),ltitle(MAXCOL),llabel(4),
     -  lsource(3),lerr_typ(2),ixprint(MAXCOL)
      data ires /0/,types1 /'c','r','l','d'/
      data label /'COMPLEX','RECEPTOR','LIGAND','DELTAS'/
      data source /'Total Energy Decomposition',
     -  'Sidechain Energy Decomposition',
     -  'Backbone Energy Decomposition'/,lsource /26,20,29/
      data err_typ /'Standard Deviation','Standard Error of the Mean'/
      data llabel /7,8,6,6/,iresrow1p /0/,iresrow2p /0/,irescol1p /0/,
     -  irescol2p /0/,lerr_typ/18,26/
c     print *,'SUMMARIZE_AMBER_CSV inpt=',inpt
      len=1000
      inpt=20
      linpstatfile=10
      call blankout(inpstatfile,1,80)
      do while (inpstatfile(linpstatfile-2:linpstatfile) .ne. 'csv')
        linpstatfile=0
        call openfile(inpt,0,'statistics',10,'old',inpstatfile,
     -    linpstatfile,notfnd,0,1,1,1,0)
        if (inpstatfile(linpstatfile-2:linpstatfile) .ne. 'csv')
     -    print *,'ERROR: file extension is not csv'
      end do
      outtabfile=inpstatfile
      outtabfile(linpstatfile-2:linpstatfile)='sum'
      iout=21
      call openfile(iout,0,' ',1,'new',outtabfile,
     -  linpstatfile,notfnd,0,1,1,1,0)
      write (iout,2000) inpstatfile(1:linpstatfile)
      ipairs=0
      ierr=0
      idecomp=0
      irepeat=0
      call blankout(line,1,len)
      line(1:1)='|'
      do while (line(1:1) .eq. '|')
        call blankout(line,1,len)
        read (inpt,1000,end=999) line
      end do
      if (line(1:9) .eq. 'idecomp =') read (line(10:11),*) idecomp
      if (idecomp .eq. 4) ipairs=1
      call blankout(atitle,1,80)
      read (inpt,1000,end=999) atitle
      call lastchar(atitle,latitle,80)
      call getskipcomma(inpt,line,len,llab,ifail)
      lsyslabel=llab
      syslabel(1:lsyslabel)=line(1:llab)
      call findname(syslabel,label,1,4,nc,lsyslabel)
      write (6,*) atitle(1:latitle)
      write (iout,*) atitle(1:latitle)
100   call quiz(ans,nt,types1(nc),'d',0,'Data type to tabulate',21,0,
     -  5,6,0)
      if (nt .ne. nc) then
        do while (line(1:llabel(nt)) .ne. label(nt)(1:llabel(nt)))
          call getskipcomma(inpt,line,len,llab,ifail)
          if (ifail .gt. 0) then
            write (6,1019) label(nt)(1:llabel(nt)),
     -        outtabfile(1:linpstatfile)
            if (irepeat .gt. 0) stop
            rewind inpt
            go to 100
          end if
        end do
      end if
      call quiz(ans,ns,'t',' ',0,'Data source to tabulate',23,0,5,6,0)
      do while (line(1:lsource(ns)) .ne. source(ns)(1:lsource(ns)))
        call getskipcomma(inpt,line,len,llab,ifail)
        if (ifail .gt. 0) then
          write (6,1019) source(ns)(1:lsource(ns)),
     -      outtabfile(1:linpstatfile)
          if (irepeat .gt. 0) stop
          rewind inpt
          go to 100
        end if
      end do
      write (iout,1017) 'type',label(nt)(1:llabel(nt))
      write (iout,1017) 'source',source(ns)(1:lsource(ns))
c     lenergy_source=llab
c     energy_source(1:lenergy_source)=line(1:llab)
      call blankout(line,1,len)
      read (inpt,1000) line
      ic=1
      ncol=0
      call findnextchar(',',line,ic,len)
      ic=ic+1
      call findnextchar(',',line,ic,len)
      ic1=ic+1
      do while (ic .lt. len)
        call findnextchar(',',line,ic,len)
c       print *,'IC=',ic
        if (ic .gt. ic1+1 .and. ic .lt. len) then
          ncol=ncol+1
          call blankout(title(ncol),1,30)
          ltitle(ncol)=ic-ic1
          title(ncol)(1:ltitle(ncol))=line(ic1:ic-1)
        end if
        ic=ic+1
        ic1=ic
        ic=ic+1
      end do
      read (inpt,*)
      ierr_mean=2
      call askyn('Do you want to tabulate the errors too',38,1,-1,
     -  ierr,0,0)
      if (ierr .gt. 0) call quiz(ans,ierr_mean,'m',' ',0,
     -  'error estimate type',19,0,5,6,000)
      write (iout,1023) err_typ(ierr_mean)(1:lerr_typ(ierr_mean))
      call askyn('Do you want to mark with m and M the extreme values',
     -  51,1,1,mark,117,4)
      if (mark .gt. 0) write (iout,2004)
      if (ipairs .eq. 0) print *,'Residue-molecule data found'
      if (ipairs .eq. 1) print *,'Residue-residue data found'
      call read_amb_data_csv(data,err,ncol,ixr1,ixr2,resnames,nresdat,
     -  ipairs,ierr_mean,inpt,MAXCOL,mxrsd,MAXREC)
      print *,'Number of properties in the input data:',ncol
      if (ipairs .eq. 0) then
        print *,'Number of residues read=',nresdat
        nres=nresdat
      else
        print *,'Number of residue pairs read=',nresdat
        nrespair=nresdat
        do ir=1,nresdat
          ixfres(ir)=0
          ixlres(ir)=0
        end do
        iresprev=0
        nres=0
        do ir=1,nresdat
          ires=ixr1(ir)
          if (ires .ne. iresprev) then
            nres=nres+1
            ixfres(nres)=ir
            if (nres .gt. 1) ixlres(nres-1)=ir-1
            iresprev=ires
          end if
        end do
        ixlres(nres)=nrespair
      end if
      if (ipairs .eq. 1) write (6,1022)
      iok=0
      do while (iok .eq. 0)
        ncolp=0
        nd=1
        do while (nd .gt. 0)
          call pickname('Index of choice (0 to finish the list)',38,
     -      title,ltitle,ncol,nd)
          if (nd .eq. 0) then
            if (ncolp .eq. 0) then
              print *,'WARNING: Nothing selected'
              call askyn('OK',2,1,-1,iok,0,0)
              if (iok .eq. 1) return
            else
              write (6,2001)
     -          (title(ixprint(ip))(1:ltitle(ixprint(ip))),ip=1,ncolp)
              call askyn('OK',2,1,+1,iok,0,0)
            end if
          else
            ncolp=ncolp+1
            ixprint(ncolp)=nd
          end if
        end do
      end do
      maxresno=ixr1(nresdat)
      call zeroiti(iresix,0,maxresno)
      do ir=1,nres
        if (ipairs .eq. 0) iresix(ixr1(ir))=ir
        if (ipairs .eq. 1) iresix(ixr2(ir))=ir
      end do
      call zeroit(colsumm,ncolp)
      do ir=1,maxresno
        minmaxlab1(ir)=' '
      end do
      if (ipairs .eq. 0) then
        if (mark .eq. 1) then
c         Establish min and max
          do ip=1,ncolp
            rmin(ip)=data(ixprint(ip),1)
            rmax(ip)=data(ixprint(ip),1)
            ixcolmin(ip)=1
            ixcolmax(ip)=1
          end do
          do ir=2,nres
            do ip=1,ncolp
              r=data(ixprint(ip),ir)
              if (r .gt. rmax(ip)) then
                rmax(ip)=r
                ixcolmax(ip)=ir
              else if (r .lt. rmin(ip)) then
                rmin(ip)=r
                ixcolmin(ip)=ir
              end if
            end do
          end do
        end if
        if (ierr .eq. 0 .or. ncolp .gt. 4)
     -    write (iout,1001)
     -     (title(ixprint(ip))(1:ltitle(ixprint(ip))),ip=1,ncolp)
        if (ierr .eq. 0) then
          do ir=1,nres
            if (mark .eq. 1) call labminmax(ir,ncolp,ixcolmin,
     -        ixcolmax,minmaxlab1,mxrsd)
            write (iout,1002) ixr1(ir),resnames(ixr1(ir))(1:3),
     -        (data(ixprint(ip),ir),minmaxlab1(ip),ip=1,ncolp)
            end do
        else
          if (ncolp .gt. 4) then
            do ir=1,nres
              if (mark .eq. 1) call labminmax(ir,ncolp,ixcolmin,
     -          ixcolmax,minmaxlab1,mxrsd)
              write (iout,1002) ixr1(ir),resnames(ixr1(ir))(1:3),
     -          (data(ixprint(ip),ir),minmaxlab1(ip),ip=1,ncolp)
              write (iout,1007)ixr1(ir),(err(ixprint(ip),ir),ip=1,ncolp)
            end do
          else
            write (iout,1004)
     -        (title(ixprint(ip))(1:ltitle(ixprint(ip))),ip=1,ncolp)
            do ir=1,nres
              if (mark .eq. 1) call labminmax(ir,ncolp,ixcolmin,
     -          ixcolmax,minmaxlab1,mxrsd)
              write (iout,1003) ixr1(ir),resnames(ixr1(ir))(1:3),
     -          (data(ixprint(ip),ir),minmaxlab1(ip),
     -           err(ixprint(ip),ir),ip=1,ncolp)
            end do
          end if
        end if
        do ir=1,nres
          do ip=1,ncolp
            colsumm(ip)=colsumm(ip)+data(ixprint(ip),ir)
          end do
        end do
        if (ierr .eq. 0) then
          write (iout,1112) (colsumm(ip),ip=1,ncolp)
        else
          write (iout,1113) (colsumm(ip),ip=1,ncolp)
        end if
      else
c       Pairs
        if (ncolp .gt. 1) then
          call getrange(ires1,ixr1(1),ires2,maxresno,incr,0,
     -      'reference residue number',24,maxresno,144)
          write (iout,1114) ires1,ires2
          call zeroit(colsum,ncolp)
c         write (6,8792) (ixprint(ip),ip=1,ncolp)
c8792     format(/,' IXPRINT=',10i4)
          do ir=1,nres
            irp=ixfres(ir)
            if (ixr1(irp) .ge. ires1 .and.  ixr1(irp) .le. ires2) then
c             Requested range found
              if (mark .eq. 1) then
c               Obtain the min/max labels for this residue
                do ip=1,ncolp
                  rmin(ip)=data(ixprint(ip),ir)
                  rmax(ip)=data(ixprint(ip),ir)
                  ixcolmin(ip)=1
                  ixcolmax(ip)=1
                end do
                if (mark .eq. 1) then
                  do irr=irp,ixlres(ir)
                    do ip=1,ncolp
                      r=data(ixprint(ip),irr)
                      if (r .gt. rmax(ip)) then
                        rmax(ip)=r
                        ixcolmax(ip)=irr-irp+1
                      else if (r .lt. rmin(ip)) then
                        rmin(ip)=r
                        ixcolmin(ip)=irr-irp+1
                      end if
                    end do
                  end do
                end if
              end if
              ir0=irp-1
              call zeroit(colsum,ncolp)
              if (ierr .eq. 1) then
                write (iout,1104)
     -            title(ixprint(1))(1:ltitle(ixprint(1))),
     -            (title(ixprint(ip))(1:ltitle(ixprint(ip))),ip=2,ncolp)
                do irr=ixfres(ir),ixlres(ir)
                  if (mark .eq. 1) call labminmax(irr-irp+1,ncolp,
     -              ixcolmin,ixcolmax,minmaxlab1,mxrsd)
                  write (iout,1107)
     -              ixr1(irr),resnames(ixr1(irr))(1:3),
     -              ixr2(irr),resnames(ixr2(irr))(1:3),
     -              (data(ixprint(ip),irr),minmaxlab1(ip),
     -              err(ixprint(ip),irr),ip=1,ncolp)
                  do ip=1,ncolp
                    colsum(ip)=colsum(ip)+data(ixprint(ip),irr)
                  end do
                end do
                write (iout,1111) ' ---- SUM',(colsum(ip),ip=1,ncolp)
                do ip=1,ncolp
                  colsumm(ip)=colsumm(ip)+colsum(ip)
                end do
              else
                write (iout,1101)
     -            (title(ixprint(ip))(1:ltitle(ixprint(ip))),ip=1,ncolp)
                do irr=ixfres(ir),ixlres(ir)
                  if (mark .eq. 1) call labminmax(irr-irp+1,ncolp,
     -              ixcolmin,ixcolmax,minmaxlab1,mxrsd)
                  write (iout,1102)
     -              ixr1(irr),resnames(ixr1(irr))(1:3),
     -              ixr2(irr),resnames(ixr2(irr))(1:3),
     -              (data(ixprint(ip),irr),minmaxlab1(ip),
     -              ip=1,ncolp)
                  do ip=1,ncolp
                    colsum(ip)=colsum(ip)+data(ixprint(ip),irr)
                  end do
                end do
                write (iout,1110) ' ---- SUM',(colsum(ip),ip=1,ncolp)
                do ip=1,ncolp
                  colsumm(ip)=colsumm(ip)+colsum(ip)
                end do
              end if
            end if
          end do
          if (ierr .eq. 0) then
            write (iout,1110) ' SUMofSUM',(colsumm(i),i=1,ncolp)
          else
            write (iout,1111) ' SUMofSUM',(colsumm(i),i=1,ncolp)
          end if
        else
c         Single property - print residue-residue matrix
          write (iout,1009) title(ixprint(1))(1:ltitle(ixprint(1)))
300       ifail12=1
          do while (ifail12 .eq. 1)
            write (6,1011) 'rows'
            call getrange(iresrow1,ixr1(1),iresrow2,maxresno,incr,0,
     -        'residue number',14,maxresno,0)
            call restoix(iresix,iresrow1,iresrow2,ix1,ix2,ifail12,
     -        maxresno)
            write (6,1011) 'columns'
            call getrange(irescol1,ixr1(1),irescol2,maxresno,incr,0,
     -        'residue number',14,maxresno,0)
            call restoix(iresix,iresrow1,iresrow2,ix1,ix2,ifail,
     -        maxresno)
            ifail12=ifail12+ifail
          end do
          if ((irescol1 .ge. iresrow1 .and. irescol1 .le. iresrow2) .or.
     -       (irescol2 .ge. iresrow1 .and. irescol2 .le. iresrow2)) then
            print *,'Residue ranges overlap'
            call askyn('Do you want to change the ranges',32,1,-1,ichng,
     -        0,0)
            if (ichng .eq. 1) go to 300
          end if
          do ir=1,nres
            if (ixr2(ir) .eq. iresrow1) iresrow1p=ir
            if (ixr2(ir) .eq. iresrow2) iresrow2p=ir
            if (ixr2(ir) .eq. irescol1) irescol1p=ir
            if (ixr2(ir) .eq. irescol2) irescol2p=ir
          end do
          if (irescol1p*irescol2p*iresrow1p*iresrow2 .eq. 0) then
            if (irescol1p .eq. 0) write (6,1015) irescol1
            if (irescol2p .eq. 0) write (6,1015) irescol2
            if (iresrow1p .eq. 0) write (6,1015) iresrow1
            if (iresrow2p .eq. 0) write (6,1015) iresrow2
            go to 300
          end if
          ntabcol=irescol2p-irescol1p+1
          ntabrow=iresrow2p-iresrow1p+1
c         write (6,*) 'irescol1,irescol2,iresrow1,iresrow2=',
c    -             irescol1,irescol2,iresrow1,iresrow2
c         write (6,*) 'irescol1p,irescol2p,iresrow1p,iresrow2p=',
c    -             irescol1p,irescol2p,iresrow1p,iresrow2p
c         write (6,*) 'nres,ntabcol,ntabrow=',nres,ntabcol,ntabrow
          if (ntabcol .gt. ntabrow) then
           print *,'Number of columns is larger than the number of rows'
            call askyn('Do you want to change the ranges',32,1,1,ichng,
     -        0,0)
            if (ichng .eq. 1) go to 300
          end if
          do ir=1,ntabcol
            ixcol(ir)=irescol1p-1+ir
          end do
          do ir=1,ntabrow
            ixrow(ir)=iresrow1p-1+ir
          end do
c          write (6,7711) 'ixr1:',(ixr1(ii),ii=1,nres)
c          write (6,7711) 'ixr2:',(ixr2(ii),ii=1,nres)
c          write (6,7711) 'ixrow:',(ixrow(ii),ii=1,ntabrow)
c          write (6,7711) 'ixcol:',(ixcol(ii),ii=1,ntabcol)
c7711      format(1x,a,/(20i4))
302       call zeroit(colsum,ntabcol)
          call zeroit(colsumm,ntabcol)
          call zeroiti(ixcolmin,0,ntabcol)
          call zeroiti(ixcolmax,0,ntabcol)
          call zeroit(tabrowsum,ntabrow)
          do irc=1,ntabcol
            colmin(irc)=999999.0
            colmax(irc)=-colmin(irc)
          end do
          do irr=1,ntabrow
            if0=ixfres(ixrow(irr))-1
            do irc=1,ntabcol
              r=data(ixprint(1),if0+ixcol(irc))
              colsum(irc)=colsum(irc)+r
              tabrowsum(irr)=tabrowsum(irr)+r
              if (r .lt. colmin(irc)) then
                colmin(irc)=r
                ixcolmin(irc)=irr
              end if
              if (r .gt. colmax(irc)) then
                colmax(irc)=r
                ixcolmax(irc)=irr
              end if
c             write (77,8972) irr,if0,ixrow(irr),irc,ixcol(irc),r
c8972         format(' irr=',i4,' if0=',i5,' ixrow(irr)=',i5,' irc=',i4,
c    -          ' ixcol(irc)=',i4,' data=',f10.6)
            end do
          end do
          write (iout,1012) (ixr2(irc),resnames(ixr2(irc))(1:3),
     -      irc=irescol1p,irescol2p),0,'SUM'
          write (iout,1013) ('-----------',irc=1,ntabcol+1)
          do irr=1,ntabrow
            if (mark .eq. 1) call labminmax(irr,ntabcol,ixcolmin,
     -        ixcolmax,minmaxlab1,mxrsd)
            if0=ixfres(ixrow(irr))-1
            write (iout,1014)
     -        ixr2(ixrow(irr)),resnames(ixr2(ixrow(irr)))(1:3),
     -        '|',(data(ixprint(1),if0+ixcol(irc)),minmaxlab1(irc),
     -        irc=1,ntabcol),tabrowsum(irr)
            if (ierr .gt. 0) write (iout,1016)
     -        '|',(err(ixprint(1),if0+ixcol(irc)),irc=1,ntabcol)
          end do
          write (iout,1013) ('-----------',irc=1,ntabcol+1)
          write (iout,1014) 0,'SUM',':',(colsum(i),' ',i=1,ntabcol)
          call askyn('Do you want to tabulate an other energy term',44,
     -      1,-1,more,0,0)
          if (more .eq. 1) then
303         call pickname('Index of choice (0 to finish the list)',38,
     -        title,ltitle,ncol,nd)
            if (nd .eq. 0) go to 303
            ixprint(1)=nd
            write (iout,*)
            write (iout,1009) title(ixprint(1))(1:ltitle(ixprint(1)))
            go to 302
          end if
        end if
      end if
      icorr=0
      call askyn(
     -  'Do you want the correlations between the properties read',56,
     -  1,-1,icorr,0,0)
      ncorrcalc=0
      if (icorr .gt. 0) then
        write (6,1020) outtabfile(1:linpstatfile),
     -    (title(i)(1:ltitle(i)),colsumm(i),i=1,ncolp)
        print *,'NOTE: correlation calculation is NOT limited to the ',
     -    'properties tabulated'
c       print *,'MAXRESSNO=',maxresno,' NRES=',nres,' IPAIRS=',ipairs
c       print *,'IXR1:'
c       write (6,8768) (ixr1(ir),ir=1,nres)
c       print *,'IXR2'
c       write (6,8768) (ixr2(ir),ir=1,nres)
c       print *,'IRESIX:'
c       write (6,8768) (iresix(ir),ir=1,maxresno)
c8768   format(i3,19i4)
        itypedone=0
        ntypepairs=0
        do while (itypedone .eq. 0)
          ntypepairs=ntypepairs+1
          call pickname('First property to correlate (zero to exit)',42,
     -      title,ltitle,ncol,ic1)
          if (ic1 .eq. 0) then
            itypedone=1
          else
            call pickname('Second property to correlate',28,title,
     -        ltitle,ncol,ic2)
            s1=0.0
            s2=0.0
            s12=0.0
            ss1=0.0
            ss2=0.0
            if (ipairs .eq. 0) then
              do i=1,nres
                s12=s12+data(ic1,i)*data(ic2,i)
                s1=s1+data(ic1,i)
                s2=s2+data(ic2,i)
                ss1=ss1+data(ic1,i)**2
                ss2=ss2+data(ic2,i)**2
              end do
              corr=(s12-s1*s2/nres)/
     -          sqrt((ss1-s1**2/nres)*(ss2-s2**2/nres))
              write (6,1005)
     -          title(ic1)(1:ltitle(ic1)),label(nc)(1:llabel(nc)),
     -          title(ic2)(1:ltitle(ic2)),label(nc)(1:llabel(nc)),
     -          corr
              write (iout,1005)
     -          title(ic1)(1:ltitle(ic1)),label(nc)(1:llabel(nc)),
     -          title(ic2)(1:ltitle(ic2)),label(nc)(1:llabel(nc)),
     -          corr
            else
              write (6,1021)
              ir1=0
              do while (ir1 .eq. 0)
                call getint('Residue number of the selected residue',38,
     -            ixr1(1),1,maxresno,irr1,00)
                ir1=iresix(irr1)
                if (ir1 .eq. 0) print *,'Non-interacting residue'
              end do
c             Find the row numbers corresponding to the resid numbers irr1,irr2
              do i=ixfres(ir1),ixlres(ir1)
                col1(i-ixfres(ir1)+1)=data(ic1,i)
                col2(i-ixfres(ir1)+1)=data(ic2,i)
              end do
              nrows=ixlres(ir1)-ixfres(ir1)+1
              do i=1,nrows
                s12=s12+col1(i)*col2(i)
                s1=s1+col1(i)
                s2=s2+col2(i)
                ss1=ss1+col1(i)**2
                ss2=ss2+col2(i)**2
              end do
              if (s12 .ne. 0.0) then
                corr=(s12-s1*s2/nrows)/
     -            sqrt((ss1-s1**2/nrows)*(ss2-s2**2/nrows))
                write (6,1010)
     -            title(ic1)(1:ltitle(ic1)),label(nc)(1:llabel(nc)),
     -            title(ic2)(1:ltitle(ic2)),label(nc)(1:llabel(nc)),
     -            irr1,resnames(irr1)(1:3),corr
                write (iout,1010)
     -            title(ic1)(1:ltitle(ic1)),label(nc)(1:llabel(nc)),
     -            title(ic2)(1:ltitle(ic2)),label(nc)(1:llabel(nc)),
     -            irr1,resnames(irr1)(1:3),corr
              else
                if (s1 .eq. 0.0) write(6,1018) title(ic1)(1:ltitle(ic1))
                if (s2 .eq. 0.0) write(6,1018) title(ic2)(1:ltitle(ic2))
              end if
            end if
          end if
        end do
      end if
      if (irepeat .eq. 0) then
        call askyn(
     -    'Do you want to tabulate a different data type or source',55,
     -    1,-1,more,0,0)
        if (more .gt. 0) then
          rewind inpt
          go to 100
        end if
      end if
      line(1:1)='*'
      do while (line(1:7) .ne. 'idecomp')
        read (inpt,1000,end=999) line
      end do
      if (line(1:9) .eq. 'idecomp =') read (line(10:11),*) idecomp
      if (idecomp .eq. 4) ipairs=1
      call blankout(atitle,1,80)
      read (inpt,1000,end=999) atitle
      call lastchar(atitle,latitle,80)
      call getskipcomma(inpt,line,len,llab,ifail)
      lsyslabel=llab
      syslabel(1:lsyslabel)=line(1:llab)
      call findname(syslabel,label,1,4,nc,lsyslabel)
      write (iout,*)
      write (iout,*) atitle(1:latitle)
      print *
      write (6,*) atitle(1:latitle)
      irepeat=1
      go to 100
999   stop
1000  format(a)
1001  format(/,'  Residue ',a10,20(1x,a9))
1002  format(i5,1x,a3,2x,20(f9.2,a1))
1003  format(i5,1x,a3,1x,4(f9.2,a1,f8.2))
1004  format(/,'  Residue ',4(a10,'    +/- '))
1005  format(' Correlation between ',a,' (',a,') and ',a,' (',a,')=',
     -  f6.4)
1007  format(i5,' +/-:',f9.2,19f10.2)
1009  format(' Residue-residue contribution table for property ',a,/)
1010  format(/,' Correlation between ',a,' (',a,') and ',a,' (',a,
     -  ') for residue ',i4,' (',a,')=',f6.4)
1011  format(' Specify residue range for the table ',a)
1012  format(13x,1000(i6,1x,a3))
1013  format(13x,1000(a))
1014  format(1x,i6,1x,a3,1x,a1,f10.2,a1,1000(f9.2,a1))
1015  format(' ERROR: residue ',i6,' is not tabulated')
1016  format(7x,'+/-: ',a1,1000(f10.2),1x)
1017  format(' Data ',a,' tabulated: ',a)
1018  format(' Sorry, can not corrwlate correlation since all ',a,
     - ' values are zero')
1019  format(1x,a,' is not found in file ',a)
1020  format(/,' Column sums of all properties in the file ',a,':',/,
     -  (1x,a20,'SUM:',f12.5))
1021  format(' Input table has residue-residue energies',/,
     -  ' Correlations will be calculated for a the different energy ',
     -  'terms',/,' of a selected residue')
1022  format(' NOTE: if you select only one property to tabulate then ',
     -  /,' a residue-residue matrix will be printed',/,' otherwise ',
     -  'separate residue-property matrices will be printed ',/,
     -  ' for each selected residue')
1023  format(' The symbols +/- stand for ',a)
1101  format(/,'  Residue - Residue ',a10,19(1x,a10))
1102  format(i5,1x,a3,i6,1x,a3,20(f10.2,a1))
1104  format(/,'  Residue - Residue ',20(a10,'    +/- '))
1107  format(i5,1x,a3,i6,1x,a3,20(f9.2,a1,f8.2))
1110  format(6x,a9,':  ',20f11.2)
1111  format(5x,a9,':  ',20(f11.2,7x))
1112  format('      SUM:',20f10.2)
1113  format('      SUM:',20(f11.2,7x))
1114  format(' Interaction tables will be generated  for residue-',
     -  'residue energies for residues in the range [',i5,',',i5,'] ',/,
     -  ' (interacting with all the rest)')
2000  format(/,' Customized tabulation of Amber residue analysis',/,
     -  ' File analyzed: ',a)
2001  format(' Column(s) to tabulate:',50(1x,a))
c2003  format(' A PDB file can be also read to provide residue names',/,
c     -  ' Residue numbers should correspond to the residue numbers ',
c     -  'read by mmpbsa.pl',/,
c     -  ' Hitting enter will skip reading/using residue names')
2004  format(' The minimum and maximum ',
     -  'values in each column are marked with m and M, resp.')
      end
