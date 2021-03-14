      subroutine explanation(ihelp,itip)
      character*60 answers,tips
      common /helplist/ init,maxans,answers(1000),
     -  linelen(1000),ifl(1000),ill(1000)
      common /tiplist/ initt,maxtips,tips(100),
     -  linelent(100),iflt(100),illt(100)
c     The answer to question iq is in lined ifl(iq)-ill(ig);
c     each line is linelen characters long
c     print *,'EXPLANATION ihelp,itip=',ihelp,itip,' init,initt=',
c    -  init,initt
      idebug=0
      if (init .eq. 0) then
        nhelp=0
c       Initialize linelen, iff, ill
        init=1
        do il=1,maxans
          call lastchar(answers(il),ilc,60)
          linelen(il)=ilc
          if (answers(il)(1:4) .eq. '****') then
            if (nhelp .gt. 0) ill(nhelp)=il-1
            if (ilc .eq. 4) go to 100
            nhelp=nhelp+1
            read (answers(il)(5:ilc),*,err=300) ih
            if (ih .ne. nhelp) go to 300
            ifl(nhelp)=il+1
          end if
        end do
100     if (idebug .gt. 0) then
          do i=1,nhelp
            print *,'i,ifl,ill=',i,ifl(i),ill(i)
            do ip=ifl(i),ill(i)
              print *,'ip,len=',ip,linelen(ip)
              print *,answers(ip)(1:linelen(ip))
            end do
          end do
        end if
        maxans=nhelp
c       print *,'NHELP=',nhelp
      end if
      if (initt .eq. 0) then
        ntips=0
c       Initialize linelent, ifft, illt
        initt=1
        do il=1,maxtips
          call lastchar(tips(il),ilc,60)
          linelent(il)=ilc
          if (tips(il)(1:4) .eq. '****') then
            if (ntips .gt. 0) illt(ntips)=il-1
            if (ilc .eq. 4) go to 200
            ntips=ntips+1
            read (tips(il)(5:ilc),*,err=300) it
            if (it .ne. ntips) go to 300
            iflt(ntips)=il+1
          end if
        end do
200     continue
c       print *,'NTIPS=',ntips
        maxtips=ntips
      end if
      if (ihelp .gt. maxans) then
        print *,'PROGRAM ERROR: invalid help number=',ihelp
        ihelp=0
      end if
      if (itip .gt. maxtips) then
        print *,'PROGRAM ERROR: invalid tip number=',itip
        itip=0
      end if
      if (ihelp .gt. 0) write (6,1000)
     -   (answers(i)(1:linelen(i)),i=ifl(ihelp),ill(ihelp))
      if (itip .gt. 0) write (6,1000)
     -   (tips(i)(1:linelent(i)),i=iflt(itip),illt(itip))
      return
300   print *,'Help data error: ',nhelp,'-th item is labeled ',ihelp
      return
1000  format(5x,a)
      end
