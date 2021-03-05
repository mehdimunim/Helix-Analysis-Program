      subroutine definelist(ansrun,nslt,nanchorr,indexn,indexo,nsegslt,
     -  segid4,iresno,molsltlim,label,llabel,iall,maxanchorlist)
      dimension indexn(nslt),indexo(nslt),iresno(nslt),
     -  molsltlim(3,nsegslt)
      character*4 segid4(nsegslt)
      character*1 ansrun
      character*(*) label
c     print *,'DEFINELIST nslt=',nslt,' NANCHORR=',nanchorr,
c    -  ' MAXANCHORLIST=',maxanchorlist,' NSEGSLT=',nsegslt
      if (ansrun .eq. 'q') then
        if (nanchorr .eq. 0) then
          print *,'No selection was made - all eligible atoms will ',
     -      'be used'
          call indexit(indexn,1,nslt,0)
          nanchorr=nslt
        end if
        return
      else if (ansrun .eq. 'a') then
c       Use all eligible solute atoms
        call indexit(indexn,1,nslt,0)
        nanchorr=nslt
        iall=1
        return
      else if (ansrun .eq. 'l') then
c       Input anchor list
        call getlist(indexo,nanchoradd,1,nslt,1,maxanchorlist)
        call trnsfi(indexn(nanchorr+1),indexo,nanchoradd)
        nanchorr=nanchorr+nanchoradd
        iall=0
      else if (ansrun .eq. 's') then
c       Input segment number
        call getint('Segment number',14,1,1,nsegslt,iseganc,0)
        nseg=molsltlim(2,iseganc)-molsltlim(1,iseganc)+1
        call indexit(indexn,nanchorr+1,nanchorr+nseg,
     -    molsltlim(1,iseganc)-nanchorr-1)
        nanchorr=nanchorr+nseg
      else if (ansrun .eq. 'x') then
c       Input anchor range
        call getrange(ifirst,1,ilast,nslt,increment,0,
     -    'anchor atoms',12,nslt,0)
        nanchoradd=ilast-ifirst+1
        call indexit(indexn,nanchorr+1,nanchorr+nanchoradd,
     -    ifirst-nanchorr-1)
        nanchorr=nanchorr+nanchoradd
        iall=0
      else if (ansrun .eq. 'n' .or. ansrun .eq. 'r') then
c       Input anchor residue list or range
        if (nsegslt .gt. 1) then
          call getint('Segment number',14,1,1,nsegslt,iseganc,0)
          write (6,2143) segid4(iseganc),(molsltlim(i,iseganc),i=1,2)
        else
          iseganc=1
        end if
        ifss=molsltlim(1,iseganc)
        ilss=molsltlim(2,iseganc)
c       Input anchor residue list or range
        if (ansrun .eq. 'n') then
c         Get list
          call getlist(indexo,nanchorres,ifss,ilss,1,maxanchorlist)
        else
c         Get range
          call getrange(ifirstres,iresno(ifss),ilastres,iresno(ilss),
     -      increment,0,'anchor residue',14,nslt,0)
          nanchorres=ilastres-ifirstres+1
          call indexit(indexo,1,nanchorres,ifirstres-1)
        end if
        do ir=1,nanchorres
          iresanc=indexo(ir)
          call findrange(iresno,ifss,ilss,iresanc,ifsr,ilsr,'residue',
     -      7,1,ifail)
          if (ifail .gt. 0) then
            write (6,2147) iresanc,ifss,ilss
          else
c           Add atoms ifsr - ilsr to the anchor list
c           print *,'IR=',ir,' ADDING ',ifsr,' - ',ilsr
            natadd=ilsr-ifsr+1
            call indexit(indexn,nanchorr+1,nanchorr+natadd,
     -        ifsr-nanchorr-1)
            nanchorr=nanchorr+natadd
c           print *,'NANCHORR=',nanchorr
          end if
        end do
c       print *,'Total number added to indexn=',nanchorr
        iall=0
      end if
      if (nanchorr .gt. maxanchorlist) then
        write (6,2144) maxanchorlist,label(1:llabel)
        nanchorr=maxanchorlist
        iall=0
        print *,'Only the first ',nanchorr,' atoms will be anchors'
        call askstop(0)
      end if
c     print *,'DEFINELIST return'
      return
2143  format(' Segment ',a,' atom range [',i6,',',i6,') selected')
2144  format(' Progam is limited to ',i5,1x,a,' anchor ',
     -  'atoms.',/,' Redimension or break up the run')
2147  format(' Residue',i6,' is not found in atom range [',i6,' - ',i6,
     -  '] - ignored')
      end
