      subroutine checkfordup(numres,n,ifres,ilres,iresno,resnames,
     -  atnames,line,index,ndupdel,ideldup,iaddnum,iaskdup,maxrec,
     -  nrescol,inamcol1,inamcol2,nnamcol,idcol)
      dimension ifres(numres),ilres(numres),iresno(n),index(n)
      character*8 resnames(numres),atnames(n),newname
      character*132 line(maxrec)
c     print *,'CHECKFORDUP NDUPDEL,IDELDUP,IADDNUM,IASKDUP=',
c    -  ndupdel,ideldup,iaddnum,iaskdup
      do ir=1,numres
        do ia=ifres(ir)+1,ilres(ir)
          do ja=ifres(ir),ia-1
            if (line(index(ia))(idcol:idcol) .eq. ' ' .and.
     -          line(index(ja))(idcol:idcol) .eq. ' ') then
              if (atnames(ia)(1:nnamcol) .eq. atnames(ja)(1:nnamcol))
     -                                                           then
                call lastchar(atnames(ia),lc,nnamcol)
                if (iaskdup .eq. 1) then
                  write (6,2000) resnames(ir)(1:nrescol),iresno(ia),
     -              ia,ja,atnames(ia)(1:nnamcol)
                  call askyn('Do you want to delete the duplicate',35,1,
     -              -1,ideldup,0,0)
                  if (ideldup .eq. 0) then
                    if (lc .lt. nnamcol) then
c                     Try adding numbers until nothing matches
                      call askyn(
     -                  'Do you want to add numbers to differentiate',
     -                  43,1,1,iaddnum,0,0)
                    end if
                  end if
                  call askyn(
     -             'Do you want to the same action with all duplicates',
     -              50,0,1,iaskdup,0,0)
                end if
                if (iaddnum .eq. 1) then
                  newname=atnames(ia)
                  nmax=10**(nnamcol-lc)-1
                  iadd=1
                  match=1
                  do while (iadd .lt. nmax .and. match .eq. 1)
                    call writeint(newname,lc+1,iadd,nnamcol-lc)
                    match=0
                    jaa=ifres(ir)
                    do while (jaa .lt. ia .and. match .eq. 0)
                      jaa=jaa+1
                      if (atnames(jaa)(1:nnamcol) .eq.
     -                    newname(1:nnamcol)) match=1
                    end do
                    iadd=iadd+1
                  end do
                  if (match .eq. 0) then
                    atnames(ia)=newname
                    line(index(ia))(inamcol1:inamcol2)=
     -                newname(1:nnamcol)
                  else
                    print *,'Sorry, could not find a different i',
     -                'new name for ',atnames(ia)(1:nnamcol)
                  end if
                else if (ideldup .eq. 1) then
                  line(index(ia))(idcol:idcol)='*'
                  ndupdel=ndupdel+1
                end if
              end if
            end if
          end do
        end do
      end do
      if (ndupdel .gt. 0)
     -  write (6,*) 'Number of duplicate names deleted=',ndupdel
      return
2000  format(' Residue ',a,i6,' atoms',2i7,' have same name:',a)
      end
