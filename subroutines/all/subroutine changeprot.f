      subroutine changeprot(resnam,resnam1,iaaconv)
c     if iaaconv=1 convert from 1-digit to 3 digit AA code
c     if iaaconv=2 convert from 3-digit to 1 digit AA code
      character*8 resnam
      character*1 aanames1,resnam1
      character*2 mmodtoamb
      character*3 aanames3
      common /atnamcon/ mmodtoamb(100),aanames1(58),aanames3(58),
     -  naanames,nnanames,nnammnames,nnames,ixwatnam
      inf=1
      inl=naanames+nnanames
      if (iaaconv .eq. 1) then
        resnam='***     '
        do i=inf,inl
          if (resnam1 .eq. aanames1(i)) then
            resnam(1:3)=aanames3(i)
            return
          end if
        end do
      else if (iaaconv .eq. 2) then
        resnam1='*'
        call leftadjustn(resnam,resnam,8)
        do i=inf,inl
          if (resnam(1:3) .eq. aanames3(i)) then
            resnam1=aanames1(i)
            return
          end if
        end do
      else
        print *,'PROGRAM ERROR: invalid iaaconv in changeprot=',iaaconv
      end if
      return
      end
