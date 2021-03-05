      subroutine quiz(char,ians,defchar,prefix,lprefix,question0,lqt0,
     -  nqincr,inpt,iout,ihelp)
      character*(*) prefix,question0
      character*60 ansline
      character*80 line,line1(100),question
      character*1 char,defchar,ans,ch(2,100)
      dimension len(100),lenans(100)
      character*60 q(100),promptlist,prompttype
      common /quizinfo/ nqst(600),lqst(600),iqfst(600),iqlst(600),
     -  maxq,init,nprompttype,promptlist(600),prompttype(600)
      common /logging/ logfile,ipredict
c     if lprefix > 0 then the question will be composed of the prefix+question0
c     if lprefix < 0 then the question will be composed of the question0+prefix
c     drop the last nqincr questions
c     print *,'QUIZ lprefix=',lprefix,' lqt0=',lqt0,' init=',init
      idebug=0
      if (init .eq. 0) then
c       Initialize nq,lq,iqfst,iqlst
        init=1
        nprompttype=0
        lqq=0
        do il=1,maxq
          call lastchar(promptlist(il),ilc,60)
          if (promptlist(il)(1:4) .eq. '****') then
            if (nprompttype .gt. 0) then
              iqlst(nprompttype)=il-1
              lqst(nprompttype)=lqq
            end if
            if (ilc .eq. 4) go to 10
            nprompttype=nprompttype+1
            prompttype(nprompttype)(1:55)=promptlist(il)(6:60)
            if (nprompttype .gt. 1) iqlst(nprompttype-1)=il-1
            iqfst(nprompttype)=il+1
          else
            if (lqq .lt. ilc) lqq=ilc
          end if
        end do
c       Quiz i has label prompptype(i) of length lqst(i)
c       List of menu items in Quiz i: promptlist(iqfst(i))-promptlist(iqlst(i))
10      if (idebug .gt. 0) then
          do i=1,nprompttype
            print *,'i,ifl,ill=',i,iqfst(i),iqlst(i)
            do ip=iqfst(i),iqlst(i)
              print *,'ip,len=',ip,lqst(i)
              print *,promptlist(ip)(1:lqst(i))
            end do
          end do
        end if
      end if
c     Find the quiz label
      iquiz=1
      do while (question0(1:lqt0) .ne. prompttype(iquiz)(1:lqt0) .and.
     -   iquiz .le. nprompttype)
        iquiz=iquiz+1
      end do
      if (iquiz .gt. nprompttype) then
        print *,'PROGRAM ERROR: quiz type ',question0(1:lqt0),
     -    ' is unknown'
        stop
      end if
      if (lprefix .eq. 0) then
        question(1:lqt0)=question0(1:lqt0)
        lqt=lqt0
      else if (lprefix .gt. 0) then
c       Combine prefix and question0
        question(1:lprefix)=prefix(1:lprefix)
        lqt=lprefix+1
        question(lqt:lqt)=' '
        question(lqt+1:lqt+lqt0)=question0(1:lqt0)
        lqt=lqt+lqt0
      else
c       Combine question0 and prefix
        question(1:lqt0)=question0(1:lqt0)
        lqt=lqt0+1
        question(lqt:lqt)=' '
        question(lqt+1:lqt-lprefix)=prefix(1:-lprefix)
        lqt=lqt-lprefix
      end if
      if (idebug .gt. 0) print *,'QUIZ nprompttype,iquiz=',
     -  nprompttype,iquiz,' nq,lq=',nq,lq
      do iq=iqfst(iquiz),iqlst(iquiz)
        q(iq-iqfst(iquiz)+1)=promptlist(iq)
      end do
      nq=iqlst(iquiz)-iqfst(iquiz)+1
      lq=lqst(iquiz)
      if (nq .gt. 40) then
        write (iout,*) '***** Redimension quiz to more than 40 items'
        stop
      end if
      nqq=nq-nqincr
      if (ihelp .gt. 0) nqq=nqq+1
      ans=' '
100   if (lqt .gt. 0) write (iout,1002) question(1:lqt)
      idebug=idebug+1
      if (idebug .eq. 4) stop
      ideffound=0
      do iq=1,nq-nqincr
c       Find signal character position (ich)
        ich=1
        do while (q(iq)(ich:ich) .ne. '<' .and. ich .lt. lq)
          ich=ich+1
        end do
        ich=ich+1
        ichgrt=ich+1
c       Find last non-blank in the menu item line
        il=lq
        do while (q(iq)(il:il) .eq. ' ' .and. il .gt. 1)
          il=il-1
        end do
        len(iq)=il
        lenans(iq)=il-2
        if (ichgrt .eq. il) lenans(iq)=lenans(iq)-1
        ch(1,iq)=q(iq)(ich:ich)
        iconvcase=1
        if (ich .lt. lq-1) then
          call findcase(q(iq)(ich+2:ich+2),icase)
          if (icase .eq. 1) iconvcase=0
        end if
        if (iconvcase .eq. 1 .and. ich .gt. 3) then
          call findcase(q(iq)(ich-2:ich-2),icase)
          if (icase .eq. 1) iconvcase=0
        end if
c       Find lower-case of signal character
        call uplow(ch(1,iq),ch(2,iq),1,noabc)
        if (noabc .eq. 1) ch(2,iq)=ch(1,iq)
c       Generate <>-less menu item
        if (ich .eq. 2) then
          line1(iq)(1:1)=ch(1,iq)
          line1(iq)(2:il-2)=q(iq)(4:il)
        else
          line1(iq)(1:ich-2)=q(iq)(1:ich-2)
          line1(iq)(ich-1:ich-1)=ch(iconvcase+1,iq)
          if (il-ich .gt. 1) line1(iq)(ich:il-2)=q(iq)(ich+2:il)
          call uplow(line1(iq)(1:1),line1(iq)(1:1),2,noabc)
        end if
c       Ask from the terminal - complete line with : and signal character
        line(1:len(iq))=q(iq)(1:len(iq))
        do ic=len(iq)+1,lq
          line(ic:ic)=' '
        end do
        do ic=len(iq)+2+mod(len(iq),2),lq-1,2
          line(ic:ic)='-'
        end do
        line(lq+1:lq+2)=': '
        line(lq+3:lq+3)=ch(2,iq)
        line(lq+4:lq+4)=ch(1,iq)
        lqw=lq+3
        if (defchar .eq. ch(2,iq)) then
          ideffound=1
          line(lq+4:lq+14)=' (default)'
          lqw=lq+14
        end if
        if (iq .lt. nqq) write (iout,1001) line(1:lqw)
        if (iq .eq. nqq) write (iout,1006) line(1:lqw)
      end do
      if (ihelp .gt. 0) then
c       Add help line
        ch(1,nqq)='?'
        ch(2,nqq)='?'
        call blankout(line,1,lq)
        line(1:4)='Help'
        do ic=6,lq-1,2
          line(ic:ic)='-'
        end do
        line(lq+1:lq+3)=': ?'
        write (iout,1006) line(1:lq+3)
      end if
c     Read and evaluate answer character
      ans=' '
      read (inpt,1000) ans
      if (defchar .ne. ' ') then
        if (ideffound .eq. 0) then
          print *,'PROGRAM ERROR: invalid default character: ',defchar
        else if (ans .eq. ' ' .or. ans .eq. '') then
          ans=defchar
        end if
      end if
      do iq=1,nq
        ians=iq
        if (ans .eq. ch(1,iq)) then
          char=ch(2,iq)
          go to 200
        else if (ans .eq. ch(2,iq)) then
          char=ans
          go to 200
        end if
      end do
      if (ihelp .gt. 0 .and. ans .eq. '?') then
        call explanation(ihelp,0)
      else
        write (iout,*) ans,': Invalid answer, pls repeat'
      end if
      go to 100
200   if (ians .ge. 1 .and. ians .le. nq) then
        lans=lenans(ians)
        write (iout,1003) line1(ians)(1:lans)
        ansline(1:lans)=line1(ians)(1:lans)
        if (logfile .gt. 0) write (logfile,1000) ans
      else if (ians .eq. nqq) then
        call explanation(ihelp,0)
        go to 100
      else
        lans=1
        ansline(1:lans)='?'
      end if
      return
1000  format(a1)
1001  format(1x,a,a)
1002  format(/,' SELECT ',a,':')
1003  format(' "',a,'" selected')
1006  format(1x,a,'  ',$)
      end
