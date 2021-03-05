      subroutine modrepats
c*****Allow for modification/addition of representative atomnames
      character*3 repats,resnam3,atnam3
      common /represent/ repats(2,50),maxrepat
      write (6,2000) ((repats(k,i),k=1,2),i=1,maxrepat)
      write (6,2001)
      call askyn('Do you want to add/modify representative atom list',
     -  50,1,-1,ians,0,0)
      if (ians .eq. 1) then
100     resnam3='   '
        call getname(resnam3,len,
     -    'Residue name to add/modify (hit ENTER to finish)',38,3,
     -    '',0,0,0,0)
        if (resnam3 .eq. '   ') then
          write (6,2000) ((repats(k,i),k=1,2),i=1,maxrepat)
          return
        end if
        atnam3='   '
        call getname(atnam3,len,
     -    'Representative atom name (3 characters maximum)',47,3,
     -    '',0,0,0,3)
        ifound=0
        do i=1,maxrepat
          if (repats(1,i) .eq. resnam3) ifound=i
        end do
        if (ifound .eq. 0) then
          maxrepat=maxrepat+1
          if (maxrepat .gt. 50) then
            print *,'Maximum number of residue names (50) exceeded'
            stop
          end if
          ifound=maxrepat
          repats(1,ifound)=resnam3
        end if
        repats(2,ifound)=atnam3
        go to 100
      end if
      return
2000  format(' Representative atoms for residues:',/,
     -  (1x,8(a3,'-',a3,2x)))
2001  format(' Representative atoms for unlisted residues: CA; ',
     -  'if not found, the 1st atom')
      end
