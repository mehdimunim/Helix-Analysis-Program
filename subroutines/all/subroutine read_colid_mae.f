      subroutine read_colid_mae(ndcol,colid,lcolid,lab,llab,line,
     -  len,iend,inpt,iout,maxcol)
      character*(*) lab,line
      character*50 colid(maxcol)
      dimension lcolid(MAXCOL)
c     print *,'READ COLID inpt,iout=',inpt,iout,' lab=',lab(1:llab)
      ndcol=1
      colid(1)='index'
      lcolid(1)=5
      ic=1
      iend=0
      do while (line(ic:ic+2) .ne. ':::')
        call blankout(line,1,200)
        read (inpt,1000,end=991,err=991) line
        call lastchar(line,lc,len)
        if (lc .gt. 1) then
          ic=1
          call nextchar(line,ic,len)
          if (line(ic:ic) .ne. '#') then
            ndcol=ndcol+1
            lcolid(ndcol)=min0(lc-ic+1,50)
c           print *,'ndcol,ic,lc=',ndcol,ic,lc
            colid(ndcol)(1:lcolid(ndcol))=line(ic:lc)
          end if
        end if
      end do
      ndcol=ndcol-1
      return
991   write (iout,2000) lab(1:llab)
      iend=1
      return
1000  format(a)
2000  format(' ERROR: run out of data while reading ',a,
     -  ' column identifiers')
      end
