      subroutine get_title_cif(inpt,linein,lineread,title,ltitle,pdbid)
      character*80 title,titles(4)
      character*200 linein
      character*4 pdbid
      dimension len(4)
c     rewind inpt
      pdbid='****'
      call zeroiti(len,0,4)
      len(4)=24
      titles(4)(1:len(4))='Title read interactively'
      ltfound=0
      do while (linein(1:10) .ne. '_atom.site' .and. ltfound .lt. 3)
        if (lineread .eq. 0) then
          call blankout(linein,1,200)
          read (inpt,1000,end=999) linein
        else
          lineread=0
        end if
        if (linein(1:5) .eq. '_data') pdbid=linein(6:9)
        if (linein(1:15) .eq. '_citation.title') then
          call blankout(linein,1,200)
          read (inpt,1000,end=999) linein
          call lastchar(linein,lc,200)
          len(1)=min0(80,lc-1)
          titles(1)(1:len(1))=linein(2:len(1)+1)
          ltfound=ltfound+1
        else if (linein(1:13) .eq. '_struct.title') then
          ic=15
          call nextstring(linein,ic,ic1,ic2,200)
          len(2)=min0(80,ic2-ic1+1)
          titles(2)(1:len(2))=linein(ic1:min0(ic2,ic1+79))
          ltfound=ltfound+1
        else if (linein(1:23) .eq. '_struct.pdbx_descriptor') then
          ic=25
          call nextstring(linein,ic,ic1,ic2,200)
          len(3)=min0(80,ic2-ic1+1)
          titles(3)(1:len(3))=linein(ic1:min0(ic2,ic1+79))
          ltfound=ltfound+1
        end if
      end do
999   write (6,2000) (i,titles(i)(1:len(i)),i=1,4)
      call getint('Choice (1-4):',12,2,1,4,in,000)
      call blankout(title,1,80)
      if (in .le. 3) then
        if (len(in) .eq. 0) then
          print *,'ERROR: this choice is empty'
          go to 999
        end if
        ltitle=len(in)
        title(1:ltitle)=titles(in)(1:ltitle)
      else
        call getname(title,ltitle,'New title',9,80,'TITLE',5,1,000,00)
      end if
1000  format(a)
2000  format(' Select the title to use:',4(/,i2,1x,a))
      end
