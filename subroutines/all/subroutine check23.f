      subroutine check23(line,index,nslt,inamcol1,inamcol2,nnamcol,
     -  maxrec)
      dimension index(nslt)
      character*132 line(maxrec)
      character*1 lcprev1,lcprev2
      character*5 atomnam,prev1,prev2
c     Look for atomnames X2,X3 without X1
      lcprev1=' '
      lcprev2=' '
      prev1(1:1)='     '
      prev2(1:1)='     '
      iasked=0
      do ia=1,nslt
        atomnam(1:nnamcol)=line(index(ia))(inamcol1:inamcol2)
        call lastchar(atomnam,lc,nnamcol)
c        write (77,2299) ia,atomnam,lc,lcprev1,lcprev2,prev1,prev2
c2299    format(i5,' anam=',a,' lc=',i1,' lcprev1,2=',a,1x,a,
c     -    ' prev1,2=',a,1x,a)
        if (atomnam(lc:lc) .eq. '3' .and.
     -      (atomnam(1:1) .eq. 'H' .or. atomnam(1:2) .eq. ' H')) then
          if ((lcprev2 .ne. '1' .or.
     -        prev2(1:lc-1) .ne. atomnam(1:lc-1)) .and.
     -        prev1(1:lc-1) .eq. atomnam(1:lc-1)) then
c           X2,X3 without X1 found
            if (iasked .eq. 0) then
              write (6,1000)
     -          atomnam(1:lc),prev1(1:lc),prev1(1:lc-1)//'1'
              call askyn('Do you want to shift H*2-H*3 to H*1-H*2',39,
     -          1,1,ishift,123,0)
              if (ishift .eq. 0) return
              iasked=1
            end if
c           Make the shift
            line(index(ia-1))(inamcol1+lc-1:inamcol1+lc-1)='1'
            line(index(ia))(inamcol1+lc-1:inamcol1+lc-1)='2'
          end if
        end if
        prev2=prev1
        prev1=atomnam
        lcprev2=lcprev1
        lcprev1=atomnam(lc:lc)
      end do
      return
1000  format(' Atom names ',a,' and ',a,' were found without ',a)
      end
