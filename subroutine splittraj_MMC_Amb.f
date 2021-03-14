      subroutine splittraj_MMC_Amb(inptrajtyp,mmctrajtyp,inpcrdtyp,
     -  nslt,naslv,c,icntrl,ifree,iconfsel,title,trtitle,indexn,indexo,
     -  indexs,lentest,maxconfsel,maxrepconf,maxrec)
      character*80 title,trtitle(32)
      dimension iconfsel(maxconfsel)
      dimension icntrl(20)
      dimension c(3,maxrec),indexn(maxrec),indexo(maxrec),
     -  indexs(maxrec),ifree(maxrec)
      real*8 xtlabc,xtlabc0
      common /boxdat/ xtlabc(6),xtlabc0(6),box(3),box0(3),edge_gen(3,3),
     -  cell0(3,27),cell(3,27),cellalt(3,27),
     -  ncell,ioppbc,noboxinfoar,noboxinfow,noboxrep,
     -  istuner,iboxtypfound,ixcrd(3),ixang(3),ixyzhex(3),
     -  ixyzhextraj(3),isizewarn
      character*200 trajnam,outfile
      common /trajname/ trajnam,outfile,ltrajnam,namleno,ifirsttraj,
     -  ifirsttraj2,ilasttraj,ilasttraj2,incrementtraj,incrementtraj2
      character*11 trajformatname
      common /trajectory/ nmmccheck,iftrajtyp(6),trajformatname(6)
      character*1 xyz
      common /axislab/ xyz(3)
c     print *,'SPLITTRAJ inptrajtyp,mmctrajtyp=',
c    -  inptrajtyp,mmctrajtyp
      if (inptrajtyp .ne. 3) then
        print *,'PROGRAM ERROR: input trajectory type is not MMC but ',
     -    trajformatname(inptrajtyp)
        return
      end if
      numsel=0
      inpt=51
      ifirst=1
      ilast=999999
      call opentraj(c,1,inpt,inptrajtyp,n,ntitltr,trtitle,
     -  inpcrdtyp,ifirst,ilast,increment,maxconf,
     -  ninconf,noutconf,natom,nfreeat,ifree,icntrl,1,mmctrajtyp,
     -  trajnam,ltrajnam,'input trajectory',16,iconfsel,numsel,
     -  0,0,0,0,icellfound,notfnd,0,0,lentest,0,ianaltyp,mx2d,
     -  maxconfsel,maxrec)
      ninconf0=ninconf
      if (noboxinfoar .eq. 0) then
        print *,'(T,P,N) run has constant number of solvents'
        return
      end if
      call askyn('Do you want to give the box size',32,0,1,nobox,0,0)
      if (nobox .eq. 0) then
        call getxyz('Box size in the ',17,' direction (A)',14,
     -    999999.0,box,1,0)
        nobox=0
        print *,'Box information is written after the first frame only'
      else
        print *,'No box information is written'
      end if
      ifail=0
      call zeroiti(indexs,0,maxrec)
      do while (ninconf .lt. ilast .and. ifail .eq. 0)
c       Read a conformation
        call readtraj(inpt,inptrajtyp,mmctrajtyp,n,naslv,nwatr,nslt,
     -    trajformatname(inptrajtyp),ntitltr,trtitle,trajnam,ltrajnam,
     -    natom,nfreeat,ifree,icntrl,c,ninconf,noutconf,
     -    increment,inpcrdtyp,ietot,etot,ifail,ifirst,ilast,iconfsel,
     -    numsel,maxrepconf,nmc,lentest,tofac,maxconf,maxconfsel,maxrec)
        numsolv=(n-nslt)/naslv
        indexs(numsolv+1)=indexs(numsolv+1)+1
      end do
      call findlim(indexs,ifst,ilst,maxrec)
      write (6,1001) ifst-1,ilst-1
      do it=ifst,ilst
c       Open output trajectory files
        if (indexs(it) .gt. 0) then
          indexo(it)=52+(it-ifst)
          lroot=ltrajnam
          if (trajnam(ltrajnam-3:ltrajnam) .eq. '.hst') lroot=lroot-4
          outfile=trajnam
          lroot=lroot+1
          outfile(lroot:lroot)='_'
          loutfile=lroot+1
          call writeint(outfile,loutfile,it-1,lenw)
          outfile(loutfile:loutfile+3)='.trj'
          loutfile=loutfile+3
          open (unit=indexo(it),status='new',form='FORMATTED',
     -      file=outfile(1:loutfile),iostat=iopen)
          if (iopen .eq. 0) then
            print *,'File ',outfile(1:loutfile),' opened'
          else
            print *,'Problem opening file ',outfile(1:loutfile)
            return
          end if
          write (title(60:79),1004) it
          write (indexo(it),1000) title
        end if
      end do
      rewind inpt
      ninconf=ninconf0
      ifail=0
      call zeroiti(indexn,ifst-1,ilst)
      do while (ninconf .lt. ilast .and. ifail .eq. 0)
c       Read a conformation
        call readtraj(inpt,inptrajtyp,mmctrajtyp,n,naslv,nwatr,nslt,
     -    trajformatname(inptrajtyp),ntitltr,trtitle,trajnam,ltrajnam,
     -    natom,nfreeat,ifree,icntrl,c,ninconf,noutconf,
     -    increment,inpcrdtyp,ietot,etot,ifail,ifirst,ilast,iconfsel,
     -    numsel,maxrepconf,nmc,lentest,tofac,maxconf,maxconfsel,maxrec)
        numsolv=(n-nslt)/naslv
        indexn(numsolv+1)=indexn(numsolv+1)+1
        iout=indexo(numsolv+1)
        write (iout,1003) ((c(k,i),k=1,3),i=1,n)
        if (indexn(numsolv+1) .eq. 1 .and. nobox .eq. 0)
     -    write (iout,1003) box0
        if (ilast .ne. 999999 .and. ilast .ge. 1000) then
          ilast10=ilast/10
          if (mod(ninconf,ilast10) .eq. 0)
     -      print *,(ninconf/ilast10)*10,'% done'
        end if
      end do
      do it=ifst,ilst
        write (6,1002) it-1,indexs(it)
        close (indexo(it))
      end do
      close (inpt)
      return
1000  format(a)
1001  format(' The number of solvents range from ',i5,' to ',i5)
1002  format(' Completed writing of trajectory of ',i5,' solvents ',
     -  'with ',i6,' configurations')
1003  format(10f8.3)
1004  format(' # of solvents=',i5)
      end
