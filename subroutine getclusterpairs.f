      subroutine getclusterpairs(npairs,iclustermem,ifstclst1,ifstclst2,
     -  ilstclst2,nslt,idebug,maxpairs,maxclustermem)
      dimension iclustermem(maxclustermem),ifstclst1(maxpairs),
     -  ifstclst2(maxpairs),ilstclst2(maxpairs)
      character*39 pairprompt
      character*80 line
      character*200 outfiletmp
c     print *,'GETCLUSTERPAIRS idebug,maxpairs,maxclustermem=',
c    -  idebug,maxpairs,maxclustermem
      listinp=55
      npairread=0
      namlenl=-1
      print *,
     -  '     Hit enter to input the cluster member list interactively'
      call openfile(listinp,1,'cluster member list',19,'old',outfiletmp,
     -  namlenl,notfnd,0,1,1,0,0)
      if (namlenl .eq. 0) then
        listinp=5
        if (npairread .eq. 0) write (6,1002)
      end if
      iread=0
      np1=-1
      do while (np1 .ne. 0)
100     if (listinp .eq. 5) then
          write (pairprompt,1000) npairread+1,1
          call getint(pairprompt,39,0,1,maxclustermem,np1,0)
        else
          read (listinp,*,err=100,end=999) np1
        end if
        if (np1 .eq. 0) go to 999
101     if (listinp .eq. 5) then
          call getname(line,lline,'Members (comma-separated):',26,80,
     -      '',0,1,0,0)
          read (line(1:lline),*,err=101,end=101)
     -     (iclustermem(iread+i),i=1,np1)
        else
          read (listinp,*,err=101,end=888)
     -      (iclustermem(iread+i),i=1,np1)
        end if
        ifstclst1(npairread+1)=iread+1
        iread=iread+np1
        ifstclst2(npairread+1)=iread+1
102     if (listinp .eq. 5) then
          write (pairprompt,1000) npairread+1,2
          call getint(pairprompt,39,1,1,maxclustermem,np2,0)
        else
          read (listinp,*,err=102,end=999) np2
        end if
103     if (listinp .eq. 5) then
          call getname(line,lline,'Members (comma-separated):',26,80,
     -      '',0,1,0,0)
          read (line(1:lline),*,err=103,end=103)
     -      (iclustermem(iread+i),i=1,np2)
        else
          read (listinp,*,err=103,end=888)
     -      (iclustermem(iread+i),i=1,np2)
        end if
        iread=iread+np2
        ilstclst2(npairread+1)=iread
        npairread=npairread+1
      end do
999   do ip=1,npairread
        if (idebug .gt. 0) then
          write (6,1003) ip,1,ifstclst2(ip)-ifstclst1(ip),
     -      (iclustermem(i),i=ifstclst1(ip),ifstclst2(ip)-1)
          write (6,1003) ip,2,ifstclst2(ip)-ifstclst2(ip)+1,
     -      (iclustermem(i),i=ifstclst2(ip),ilstclst2(ip))
        end if
        do i=ifstclst1(ip),ifstclst2(ip)-1
          if (iclustermem(i) .lt. 1 .or.  iclustermem(i) .gt. nslt) then
            write (6,1004) ip,1,iclustermem(i),nslt
          end if
        end do
        do i=ifstclst2(ip),ilstclst2(ip)
          if (iclustermem(i) .lt. 1 .or.  iclustermem(i) .gt. nslt) then
            write (6,1004) ip,2,iclustermem(i),nslt
          end if
        end do
      end do
      npairs=npairread
      return
888   print *,'ERROR: invalid list for cluster',ic12,
     -  ' pair #',npairread+1
      return
1000  format('Pair ',i3,' number of members in cluster',i2)
1002  format(' Type 0 to finish the list')
1003  format(' Pair',i3,' cluster',i2,':',(10i5))
1004  format(' ERROR: pair',i3,', cluster',i2,' member (',i6,
     -  ') is outside the',' range [1,',i6,']')
      end
