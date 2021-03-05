      subroutine opentraj(c,irept,inpt,inptrajtyp,n,ntitltr,trtitle,
     -  inpcrdtyp,ifirst,ilast,increment,maxconf,ninconf,noutconf,
     -  natom,nfreeat,ifree,icntrl,iaskalways,mmctrajtyp,inpfile,
     -  namlent,prompt,lprompt,iconfsel,numsel,iverb,notitprint,
     -  icont_traj,idatapath,icellfound,notfnd,iconfirmname,
     -  iopen_rep,lentest,iout_lim,ianaltyp,mx2d,maxconfsel,maxrec)
      dimension c(3,maxrec),ifree(maxrec)
      character*(*) prompt
      character*80 trtitle(32)
      character*200 inpfile
      dimension iconfsel(maxconfsel)
      character*4 chhd
      character*11 trajformatname
      common /trajectory/ nmmccheck,iftrajtyp(6),trajformatname(6)
      real*8 xtlabc,xtlabc0
      common /boxdat/ xtlabc(6),xtlabc0(6),box(3),box0(3),edge_gen(3,3),
     -  cell0(3,27),cell(3,27),cellalt(3,27),
     -  ncell,ioppbc,noboxinfoar,noboxinfow,noboxrep,
     -  istuner,iboxtypfound,ixcrd(3),ixang(3),ixyzhex(3),
     -  ixyzhextraj(3),isizewarn
      real*8 xtlabct
      dimension xtlabct(7),z(10),zz(10),edg(3),icntrl(20)
      character*5 crdext
      common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
     -  iommod,iommc,iommc4,iogro,iomol2,iomae,iocif,ioxxx,ioins,
     -  ionxyz,iosxyz,iosxyzrq,iograsp,iofull,lext(19),crdext(19)
      character*12 recordn
c     print *,'OPENTRAJ iverb,notitprint,irept,iopen_rep=',
c    -  iverb,notitprint,irept,iopen_rep
c     print *,'OPENTRAJ inptrajtyp,mmctrajtyp=',inptrajtyp,mmctrajtyp,
c    -  ' ifirst,ilast=',ifirst,ilast
c     print *,'OPENTRAJ icont_traj=',icont_traj,' IOUT_LIM=',iout_lim
      call blankout(recordn,1,12)
      icellfound=0
      icrorg=0
      if (inptrajtyp .eq. 0) then
        print *,'PROGRAM ERROR: missing trajectory type for ',
     -    prompt(1:lprompt)
        stop
      end if
      ntry=0
      iform=iftrajtyp(inptrajtyp)
      if (inptrajtyp .eq. 3 .and. mmctrajtyp .eq. 1) iform=2
      noecho=0
      if (iverb .eq. 0) noecho=1
      if (iconfirmname .gt. 0 .and. namlent .gt. 0) then
        print *,'Input trajectory file name=',inpfile(1:namlent)
        call askyn('Do you want to change it',24,1,-1,newn,0,0)
        if (newn .gt. 0) namlent=0
        iconfirmname=0
      end if
100   call openfile(inpt,0,prompt,lprompt,'old',inpfile,namlent,notfnd,
     -  idatapath,iform,1,noecho,0)
      if (notfnd .eq. 1) return
      rewind inpt
      if (icont_traj .le. 1) then
        ninconf=0
        noutconf=0
      end if
      ntitltr=0
      if (inptrajtyp .eq. 1) then
c       Open a Charmm trajectory
        recordn='chhd,icntrl '
        read (inpt,err=1098,end=1098) chhd,icntrl
        if (chhd .ne. 'CORD') then
          print *,'ERROR:header label is not CORD: ',chhd
          go to 1099
        end if
        recordn='title       '
        read (inpt,err=1098,end=1099) ntitltr
c       print *,'Number of title lines=',ntitltr
        if (ntitltr .lt. 1 .or. ntitltr .gt. 32) then
          print *,'WARNING: invalid number of title lines:',
     -      ntitltr
          ntitltr=0
        else
          do i=1,ntitltr
            call blankout(trtitle(i),1,ntitltr)
          end do
          rewind inpt
          read (inpt)
          read (inpt,err=1098,end=1099) ntitltr,
     -      (trtitle(i),i=1,ntitltr)
        end if
        if (notitprint .eq. 0) then
          if (trtitle(1)(1:21) .eq. 'Created by DCD plugin') then
            print *,'Trajectory written by VMD'
          else
            if (ntitltr .gt. 1) then
              if (trtitle(ntitltr-1)(1:8) .eq. 'Simulaid') then
                print *,trtitle(ntitltr-1)
                print *,trtitle(ntitltr)
              end if
            end if
            if (ntitltr .eq. 2) then
c             Check for NAMD origin
              ic=0
              do while (ic .lt. 65 .and.
     -                  trtitle(1)(ic+1:ic+15) .ne. 'CREATED BY NAMD')
                ic=ic+1
              end do
              if (trtitle(1)(ic+1:ic+15) .eq. 'CREATED BY NAMD')
     -          print *,'Trajectory written by NAMD'
            end if
          end if
        end if
        read(inpt,err=1098,end=1099) natom
        if (iverb .gt. 0 .and. natom .ne. n) then
          write (6,1020) natom,n
          call askstop(1)
        end if
        nfreeat = natom - icntrl(9)
        if (chhd .ne. 'CORD') then
          write (6,1003) chhd
          namlent=0
          go to 100
        end if
        if (iverb .gt. 0)
     -    write (6,1001) icntrl(1),nfreeat,icntrl(5),icntrl(20)
        if (n .ne. natom) then
          write (6,1002) n,natom
          call askstop(0)
          n=natom
        end if
c       Free atom array
        recordn='ifree       '
        if (icntrl(9) .gt. 0)
     -    read (inpt,err=1098,end=1099) (ifree(i),i=1,nfreeat)
        recordn='cell info1'
        if (lentest .eq. 1) then
c         Test for presence of cell information
          icell=0
          xtlabct(7)=999999.9
          read (inpt,err=120,end=120) xtlabct
120       if (xtlabct(7) .eq. 999999.9) icell=1
        else
          icell=icntrl(11)
          if (icell .eq. 1) icell=2
        end if
c       Rewind and skip back (to be on the safe side)
        rewind inpt
        read(inpt)
        read(inpt)
        read(inpt)
        if (icntrl(9) .gt. 0) read(inpt)
        if (icell .gt. 0) then
c         Cell information was found
          read (inpt,err=1098,end=1098) (xtlabc(i),i=1,6)
          noboxinfoar=0
          icellfound=1
          iboxtypfound=1
          call trnsfrd(xtlabc0,xtlabc,6)
          if (iverb .gt. 0 .or. ioppbc .eq. 4 .or. ioppbc .eq. 5) then
            write (6,1007) 'size',(xtlabc0(ixcrd(k)),k=1,3)
            write (6,1007) 'angle/cos',(xtlabc0(ixang(k)),k=1,3)
          end if
          if (ioppbc .eq. 1 .or. ioppbc .eq. 2) then
            do k=1,3
              edg(k)=xtlabc0(ixcrd(k))
            end do
            icrorg=1
          else if (ioppbc .eq. 5) then
c           Skewed hexagon
            edg(1)=xtlabc0(ixcrd(3))
            edg(2)=xtlabc0(ixcrd(2))
            edg(3)=xtlabc0(ixcrd(1))
            icrorg=1
          end if
          if (lentest .gt. 0) then
c           Test for cell data before 2nd config
            recordn='1st config'
            read(inpt,err=1098,end=1099)
            read(inpt,err=1098,end=1099)
            read(inpt,err=1098,end=1099)
            xtlabct(7)=999999.9
            recordn='cell info2'
            read(inpt,err=130,end=130) xtlabct
130         if (xtlabct(7) .eq. 999999.9) icell=2
          end if
          if (icell .eq. 2) then
            if (iverb .gt. 0)
     -        print *,'Cell information is read for each configuration'
            iboxtypfound=2
          else
            if (iverb .gt. 0)
     -        print *,'Cell information is read only at start'
          end if
          if (icntrl(11) .eq. 0) then
            print *,'Suspect: icntrl(11)=',icntrl(11)
          end if
        else
          if (iverb .gt. 0) then
            print *,'No cell information'
            if (icntrl(11) .ne. 0)
     -        print *,'Suspect: icntrl(11)=',icntrl(11)
          end if
          noboxinfoar=1
        end if
        icntrl(11)=icell
c       Reposition trajectory
        rewind inpt
        read(inpt)
        read(inpt)
        read(inpt)
        if (icntrl(9) .gt. 0) read(inpt)
      else if (inptrajtyp .eq. 2) then
c       Open an Amber trajectory
        if (nmmccheck .eq. 0 .and. inpcrdtyp .eq. iommc) then
          write (6,1021) n
          call askyn('Do you want to change it (to include solvents)',
     -      46,1,-1,newn,0,0)
          if (newn .eq. 1) call getint('New number of atoms',19,1,1,
     -      maxrec,n,0)
          nmmccheck=1
        else if (iverb .gt. 0) then
          write (6,1021) n
        end if
        recordn='title       '
        call blankout(trtitle(1),1,80)
        read (inpt,1000,err=1098,end=1099) trtitle(1)
        ntitltr=1
        if (iverb .gt. 0) write (6,1010) trtitle(1)
c       Find first box info record
        recordn='first struc.'
        nleft=3*n
        do while (nleft .gt. 0)
          nread=min0(10,nleft)
          read (inpt,1012,err=1098,end=1097) (z(i),i=1,nread)
          nleft=nleft-10
        end do
        noboxinfoar=1
        iboxtypfound=1
        call zeroit(z,10)
        read (inpt,1012,err=1098,end=1097) (z(i),i=1,10)
        call countzeros(z,10,nzero1)
c       Check if box info was present after the 2nd config.
        nzero2=0
        if (nzero1 .ge. 7) then
          icellfound=1
          noboxinfoar=0
c         Save the first box and skip to the end of the second frame
          call trnsfr(box0,z,3)
c         See if there is second box info
          recordn='secnd struc.'
          nleft=3*n
          do while (nleft .gt. 0)
            nread=min0(10,nleft)
            read (inpt,1012,err=1098,end=1097) (zz(i),i=1,nread)
            nleft=nleft-10
          end do
          call zeroit(zz,10)
          read (inpt,1012,err=1098,end=1097) (zz(i),i=1,10)
          call countzeros(zz,10,nzero2)
          if (nzero2 .ge. 7) then
            noboxrep=0
            if (iverb .gt. 0) print *,'Box information was found after',
     -        ' each structure'
          else if (iverb .gt. 0) then
            print *,'NOTE: Box information appears to be missing after',
     -        ' the first structure'
            iboxtypfound=2
          end if
        else if (iverb .gt. 0) then
          print *,'No box information was found at all'
          iboxtypfound=0
        end if
        if (nzero1 .gt. 7 .and. iverb .gt. 0)
     -    write (6,1006) 'First',(z(k),k=1,3)
        if (nzero2 .gt. 7 .and. iverb .gt. 0)
     -    write (6,1006) 'Second',(zz(k),k=1,3)
c       Reposition trajectory
        rewind inpt
        read (inpt,1000) trtitle(1)
      else if (inptrajtyp .ge. 3 .and. inptrajtyp .le. 5) then
        ntitltr=0
        if (inptrajtyp .eq. 3) then
c         MMC trajectory
          if (noboxinfoar .eq. -1 .and. mmctrajtyp .ne. 1)
     -      call askyn('Is the trajectory from a (T,P,N) run',
     -        36,0,-1,noboxinfoar,0,0)
          if (mmctrajtyp .eq. 1) then
            call binhst_type(inpt,c,ibox,istuner,ieof,6,maxrec)
            if (ieof .gt. 0) call askstop(1)
            if (iopen_rep .eq. 0) then
              if (ibox .gt. 0) print *,'Box information was found'
              if (istuner .gt. 0) print *,'Tuning information was found'
              if (ibox+istuner .eq. 0)
     -          print *,'Neither tuning nor box information was found'
            end if
            if (noboxinfoar .ne. -1 .and. noboxinfoar .ne. 1-ibox) then
              write (6,1008)
              call askstop(1)
            end if
            noboxinfoar=1-ibox
            read (inpt,err=1098,end=1999) nwat,natoms
            read (inpt,err=1098,end=1999)
            read (inpt,err=1098,end=1999) box0
1999        rewind inpt
            if (iverb .gt. 0) write (6,1011) natoms,nwat
            if (iverb .gt. 0 .and. natoms .ne. n) then
              write (6,1020) natoms,n
              call askstop(1)
            end if
          end if
          icellfound=1-noboxinfoar
          if (noboxinfoar .eq. 0) iboxtypfound=2
        end if
        call blankout(trtitle(1),1,80)
        trtitle(1)(1:27)='Trajectory generated by MMC'
      else if (inptrajtyp .eq. 6) then
c       Amber CDF
        print *,'Under implementation'
        stop
      else if (inptrajtyp .gt. 6) then
        print *,'Trajectory input type ',inptrajtyp,' is not yet done'
        stop
      end if
      if (iverb .gt. 0) then
        write (6,1013) trajformatname(inptrajtyp),inpfile(1:namlent)
        if (ntitltr .gt. 0)
     -    write (6,1014) (trtitle(i)(1:79),i=1,ntitltr)
      end if
      if (irept .eq. 1) then
        if (icrorg .eq. 1) then
          if (ioppbc .lt. 0) then
c           Set PBC to cubic, ask for size
            print *,'PBC type is unknown - it is set to cubic'
            ioppbc=1
            npbc=1
            call pbcsize(ioppbc,edg,npbc)
          end if
c         Recreate the cell from the first frame's cell size
          call crorgn(edg,edge_gen,ioppbc,3,ncell,cell,cellalt,
     -      ixyzhex,rinscr,rcirc)
          call trnsfr(cell0,cell,3*ncell)
          write (6,1016) edg
        else if (ncell .gt. 0) then
          call trnsfr(cell0,cell,3*ncell)
          if (isizewarn .eq. 1) write (6,1015)
          isizewarn=0
        end if
      end if
      if (icont_traj .eq. 0) then
        if (irept .eq. 1 .or. iaskalways .eq. 1) then
          if (numsel .eq. 0) then
c           See if a list is given for snapshots to be read
            call askyn('Do you have a list of configurations to read',
     -        44,1,-1,igetlist,0,0)
            if (igetlist .eq. 1) then
              call getlist(iconfsel,numsel,1,999999,1,maxconfsel)
              maxconf=numsel
            end if
          else
            maxconf=numsel
          end if
          if (numsel .eq. 0) then
            idefl=ilast
            ideff=max0(1,ifirst)
            if (inptrajtyp .eq. 1) idefl=icntrl(1)
200         call getrange(ifirst,ideff,ilast,idefl,increment,1,
     -        'structure to use from trajectory',32,0,0)
            if (ilast .eq. 0) ilast=999999
            maxconf=(ilast-ifirst)/increment+1
            if ((ianaltyp .eq. 23 .or. ianaltyp .eq. 24).and.
     -           maxconf .gt. mx2d) then
              write (6,1022) maxconf,mx2d,mx2d
              go to 200
            end if
            print *,'Number of configurations to use=',maxconf
          end if
          if (iout_lim .gt. 0) then
            if (numsel .gt. 0) then
              write (iout_lim,1009) (iconfsel(i),i=1,numsel)
            else if (ianaltyp .ne. 23 .and. ianaltyp .ne. 24) then
              write (iout_lim,1017) ifirst,ilast,increment
            end if
          end if
        end if
      end if
      icntrl(3)=icntrl(3)*increment
      return
1097  print *,'Amber trajectory must contain at least 3 structures'
      return
1098  write (6,1004) trajformatname(inptrajtyp),'invalid'
      print *,'TRTITLE(1)=',trtitle(1)
      print *,'Problem occured while trying to read ',recordn
      if (inptrajtyp .eq. 1) then
        print *,'chhd=',chhd
        print *,'icntrl=',icntrl
        print *,'ntitltr=',ntitltr
      end if
      namlent=0
      ntry=ntry+1
      if (ntry .lt. 3) go to 100
      stop
1099  write (6,1004) trajformatname(inptrajtyp),'incomplete'
      if (recordn(1:2) .ne. '  ')
     -  print *,'Problem occured while trying to read ',recordn
      namlent=0
      if (inptrajtyp .eq. 1) write (6,1005)
      write (6,1019) inpfile(1:namlent)
      ntry=ntry+1
      if (ntry .lt. 3) go to 100
      stop
1000  format(a)
1001  format(' Number of data sets: ',i7,/,
     -  ' Number of free atoms=',i7,' Number of fixed atoms=',i7,/,
     -  ' Charmm version=',i4)
1002  format(' WARNING: The number of atoms inputted (',i8,
     -  ') is different from',/,10x,'the number of atoms on the ',
     -  'Charmm trajectory file (',i8,')',/,10x,'The trajectory ',
     -  'file value will be used')
1003  format(' WARNING: Charm trajectory file 4-character ',
     -  'header is incorrect: ',a4,' should be CORD',/,
     -  ' Make sure that Simulaid is compiled with the same byte order',
     -  ' as the DCD file')
1004  format(' ERROR: ',a,' trajectory file entry is ',a)
1005  format(' Make sure Simulaid was compiled with the same byte ',
     -  'order',/,' as the trajectory file was run')
1006  format(' WARNING: ',a,' box information includes zero edges:',
     -  3f7.2)
1007  format(' Initial PBC box ',a,' information:',3f12.6)
1008  format(' Conflicting box information - probable progam error')
1009  format(' Frame numbers selected to be used:',/,(10i6))
1010  format(' Amber trajectory file title:',/,1x,a)
1011  format(' MMC binary trajectory - number of atoms=',i10,
     -  ' Number of solvents=',i8)
1012  format(10f8.3)
1013  format(1x,a,' trajectory file ',a,' opened')
1014  format(' Title:',/,(1x,a))
1015  format(' Make sure that the PBC size read is IDENTICAL to ',
     -  'the size of the first frame')
1016  format(' Cell edges are set from the trajectory to',/,3f10.5)
1017  format(' First frame to use=',i6,/,' Last frame to use=',i6,/,
     -  ' Increment=',i5)
1019  format(' Trajectory file opened: ',a)
1020  format(' WARNING: number of atoms found (',i8,') differs from ',
     -  'the',/,10x,'number of atoms expected (',i8,')')
1021  format(' Opening Amber trajectory, number of atoms per structure',
     -  '=',i8)
1022  format(' ERROR: Number of frames to analyze (',i6,') exceeds ',
     -  'the limit (',i5,')',/,' Increase the inrement or recompile',
     -  ' with the parameter MAX2D > ',i5,/)
      end
