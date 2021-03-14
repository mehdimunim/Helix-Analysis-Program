      subroutine openfile(iunit,iswitch,prompt,lenprompt,mode,
     -  filename,namlen,notfound,idatapath,ib,nosys,noecho,ioverall)
      character*(*) prompt,filename
      character*3 mode
      character*100 datapath,datapaths
      common /environment/ npaths,ldatapath,ldatapaths(5),
     -  datapath,datapaths(5)
      common /logging/ logfile,ipredict
      character*11 formnam
      character*80 promptq
      dimension lenform(2),formnam(2)
      data formnam /'formatted  ','unformatted'/,lenform /9,11/
c     iswitch=1: Switch to terminal input if namlen=0
c     ib: 1/2: formatted/unformatted
c     nosys=1: don't try to look into datapath
c     ioverall=1: no matter what, overwrite without asking existing file
c     idatapath=1: when an 'old' file is not found, try it in the datapath dir.
c     idatapath=2: when an 'old' file is not found, just returnt nofund=1
c     idatapath=3: when an 'old' file is not found, ask for a new file name
c     print *,'OPENFILE ipredict=',ipredict,' fn=',filename(1:namlen)
c     print *,'OPENFILE ib,f,=',ib,formnam(ib),' fn=',filename(1:namlen)
c     print *,'OPENFILE namlen=',namlen
      notfound=0
      iover=0
100   if (namlen .le. 0) then
        write (promptq,1000) prompt(1:lenprompt)
        lprompt=lenprompt+17
        call getname(filename,namlen,promptq,lprompt,200,'',0,0,0,0)
c       print *,'Name read:',filename(1:namlen)
        if (namlen .eq.  0) then
          if (iswitch .eq. 1) then
            print *,'No name was provided - switching to terminal input'
            return
          else
            go to 100
          end if
        end if
      end if
      if (ipredict .eq. 1 .and. mode .eq. 'new' .and.
     -    ioverall .eq. 0) then
        print *,'Opening file ',filename(1:namlen)
        call askyn('If the file exists, do you want to overwrite it',47,
     -    1,-1,iover,0,0)
      end if
      lenfrm=lenform(ib)
110   if (mode .eq. 'new') then
        open(unit=iunit,status='new',file=filename(1:namlen),
     -    iostat=iopen,form=formnam(ib)(1:lenfrm))
      else
        open(unit=iunit,status='old',file=filename(1:namlen),
     -    iostat=iopen,form=formnam(ib)(1:lenfrm))
      end if
      if (iopen .ne. 0) then
        if ((iover+ioverall) .gt. 0 .and. mode .eq. 'new') then
          open(unit=iunit,status='old',file=filename(1:namlen),
     -      iostat=iopen,form=formnam(ib)(1:lenfrm))
          if (iopen .gt. 0) then
            write (6,1004) 'overwriting',filename(1:namlen)
            stop
          end if
          return
        end if
        if (idatapath .eq. 2) then
c         Just return notfound=1
          notfound=1
          return
        else if (idatapath .eq. 3 .and. mode .eq. 'old') then
          write (6,1004) 'opening',filename(1:namlen)
          namlen=0
          go to 100
        end if
        if (mode .eq. 'old') then
          if (idatapath .eq. 1 .and. nosys .eq. 0 .and.
     -        ldatapath .gt. 0) then
c           Try datapath
            filename=datapath(1:ldatapath)//'/'//filename(1:namlen)
            namlen=namlen+ldatapath+1
            open(unit=iunit,status='old',file=filename(1:namlen),
     -        iostat=iopen,form=formnam(ib)(1:lenfrm))
            if (iopen .eq. 0) then
              filename(1:namlen)=filename(ldatapath+2:ldatapath+namlen)
              go to 200
            end if
          end if
          write (6,1004) 'opening',filename(1:namlen)
          namlen=-1
          go to 100
        else
          if (ipredict .eq. 0) then
            write (6,1004) 'opening',filename(1:namlen)
            call askyn(
     -        'Do you want to overwrite it if it already exists',48,
     -        1,-1,iover,0,0)
          end if
          if (iover .eq. 1) then
            print *,'Overwriting file ',filename(1:namlen)
            open(unit=iunit,status='old',file=filename(1:namlen),
     -        iostat=iopen,form=formnam(ib)(1:lenfrm))
            if (iopen .eq. 0) then
              close (iunit,status='delete')
              go to 110
            end if
            write (6,1002) filename(1:namlen)
          end if
          namlen=-1
          go to 100
        end if
      end if
200   if (namlen .eq. 0) namlen=lenchar(filename,1,200)
      if (noecho .eq. 0) write (6,1001) filename(1:namlen),
     -  formnam(ib)(1:lenfrm),iunit
      return
1000  format('Name of the ',a,' file')
1001  format(' File ',a,' (',a,') opened on unit ',i3)
1002  format(' There is still a problem opening file ',a)
1004  format(' Problem ',a,' file ',a)
      end
