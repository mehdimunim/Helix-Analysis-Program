      subroutine read_rrdist(rr,inpt,lab,llab,irefres1,irefres2,
     -  inegres1,inegres2,inpmat,linpmat,ifail,mxdim)
      dimension rr(mxdim,mxdim)
      character*(*) lab
      character*80 inpmat,prompt,line
      ifail=0
      prompt(1:llab)=lab(1:llab)
      prompt(llab+1:llab+24)=' residue distance matrix'
      lprompt=llab+24
      linpmat=0
      call openfile(inpt,0,prompt,lprompt,'old',inpmat,linpmat,notfnd,
     -  0,1,1,0,0)
      if (inpmat(linpmat-3:linpmat) .ne. '.rsd' .and.
     -    inpmat(linpmat-3:linpmat) .ne. '.rsm') write (6,1002)
      istart=1
      nresfound=0
      do while (.true.)
        read (inpt,1000,end=999) line
        if (line(1:6) .eq. ' Res #') then
          nresfound=nresfound+1
          read (line(07:12),*,end=888,err=888) irefres
          read (line(26:31),*,end=888,err=888) inegres
          if (istart .eq. 1) then
            istart=0
            irefres1=irefres
            inegres1=inegres
          end if
          read (line(43:50),*,end=888,err=888)
     -      rr(irefres-irefres1+1,inegres-inegres1+1)
        end if
      end do
999   close (inpt)
      if (nresfound .gt. 0) then
        irefres2=irefres
        inegres2=inegres
        write (6,1001) irefres1,irefres2,inegres1,inegres2,nresfound
      else
        write (6,*) 'ERROR: No residue distance record was found'
        ifail=1
      end if
      return
888   print *,'ERROR reading record'
      print *,line(1:66)
      close (inpt)
      ifail=1
      return
cRes #     1(ILE ) - Res #     1(ILE ): <d>=  0.0000 A sd=  0.0000
1000  format(a)
1001  format(' Residue ranges found: [',i6,',',i6,'] and [',i6,',',
     -  i6,']',/,' Number of distance records read=',i4)
1002  format(' WARNING: input file extension is not .rsd or .rsm')
      end
