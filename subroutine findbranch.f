      subroutine findbranch(iachir,inchir,i3chir,c,nbranch,ia_branch,
     -  ia_predef,r_predef,ba_predef,ta_predef,nneig,nneig_sav,ineig,
     -  iused,iparent,loopmem,nslt,atnames,latnam,resnames,
     -  lresnam,ixres,iresno,ifail,maxneig,maxrsd,maxat,maxbranch)
      dimension c(3,maxat),ia_branch(maxbranch),ia_predef(3,maxbranch),
     -  r_predef(maxbranch),ba_predef(maxbranch),ta_predef(maxbranch),
     -  nneig(maxat),nneig_sav(maxat),ineig(maxneig,maxat),iused(maxat),
     -  iparent(maxat),loopmem(maxbranch),ixres(maxat),iresno(maxat)
      character*8 atnames(maxat),resnames(maxrsd)
c     Atom ia_branch(i) has predecessors ia_predef(*,i),with distance,
c     angle and torsion r_predef(i),ba_predef(i),ta_predef(i)
      ifail=0
      nbranch=0
      if (nneig(inchir) .eq. 1) return
100   call zeroiti(iused,0,nslt)
      iused(iachir)=1
      iused(inchir)=1
c     move iachir to be the last neighbor of inchir
      in0=0
      do in=1,nneig(inchir)
        if (ineig(in,inchir) .eq. iachir) in0=in
      end do
      if (in0 .eq. 0) then
        print *,'PROGRAM ERROR in findbranch: iachir is not among the ',
     -    'neighbors of inchir'
        ifail=1
        return
      end if
      iparent(inchir)=iachir
      iparent(iachir)=i3chir
      nbranch=1
      ia_branch(nbranch)=inchir
      iused(ia_branch(nbranch))=1
      itry=1
      do while (itry .le. nbranch)
        iatry=ia_branch(itry)
        do while(nneig(iatry) .gt. 0)
          ianext=ineig(nneig(iatry),iatry)
          if (ianext .eq. iachir .and. iatry .ne. inchir) then
            print *,'Branch on atom ',inchir,' is a ring'
c           write (6,2001) ianext,atnames(ianext)(1:latnam),
c    -        resnames(ixres(ianext))(1:lresnam),iresno(ianext)
c           write (6,2001) iatry,atnames(iatry)(1:latnam),
c    -        resnames(ixres(iatry))(1:lresnam),iresno(iatry)
            loopmem(1)=ianext
            loopmem(2)=iatry
            nmem=2
            ip=iparent(iatry)
            do while (ip .ne. inchir)
c             write (6,2001) ip,atnames(ip)(1:latnam),
c    -          resnames(ixres(ip))(1:lresnam),iresno(ip)
              nmem=nmem+1
              loopmem(nmem)=ip
              ip=iparent(ip)
            end do
c           write (6,2001) inchir,atnames(inchir)(1:latnam),
c    -        resnames(ixres(inchir))(1:lresnam),iresno(inchir)
            nmem=nmem+1
            loopmem(nmem)=inchir
c           Open ring, restart
            print *,'Have to break one bond in the ring'
            write (6,2002) (i,loopmem(i),atnames(loopmem(i))(1:latnam),
     -        resnames(ixres(loopmem(i)))(1:lresnam),
     -        iresno(loopmem(i)),loopmem(i+1),
     -        atnames(loopmem(i+1))(1:latnam),
     -        resnames(ixres(loopmem(i+1)))(1:lresnam),
     -        iresno(loopmem(i+1)),i=1,nmem-1)
            call getint('Bond # to break (1/2/3/...)',27,nmem-1,1,
     -        nmem-1,ibreak,000)
            i1=loopmem(ibreak)
            i2=loopmem(ibreak+1)
            write (6,2003) i1,atnames(i1)(1:latnam),
     -        resnames(ixres(i1))(1:lresnam),iresno(i1),
     =        i2,atnames(i2)(1:latnam),
     -        resnames(ixres(i2))(1:lresnam),iresno(i2)
            print *,'NOTE: swapped ring will be distorted!'
            call trnsfi(nneig,nneig_sav,nslt)
            call breakbond0(i1,i2,nslt,nneig,ineig,ifail,maxneig)
            go to 100
          end if
          if (iused(ianext) .eq. 0) then
            iused(ianext)=1
            nbranch=nbranch+1
            if (nbranch .gt. maxbranch) then
              write (6,2000) maxbranch
              call trnsfi(nneig,nneig_sav,nslt)
              ifail=1
              return
            end if
            ia_branch(nbranch)=ianext
            iparent(ianext)=iatry
c           print *,'Atom ',ianext,' added to ia_branch; nbr=',nbranch
          end if
          nneig(iatry)=nneig(iatry)-1
        end do
        itry=itry+1
      end do
c     List is complete - set its properties
      do ia=2,nbranch
c       Don't use the first atom on the list - it is swapped separately
        iaa=ia-1
        ian=ia_branch(ia)
        ia_branch(iaa)=ian
        ia1=iparent(ian)
        ia2=iparent(ia1)
        ia3=iparent(ia2)
c       print *,'IAN,IA1,2,3=',ian,ia1,ia2,ia3
        ia_predef(1,iaa)=ia1
        ia_predef(2,iaa)=ia2
        ia_predef(3,iaa)=ia3
        r_predef(iaa)=sqrt(dist2(c(1,ian),c(1,ia1)))
        ba_predef(iaa)=angleijk(c,maxat,ian,ia1,ia2,6)
        ta_predef(iaa)=dihangl(c,ia3,ia2,ia1,ian,0,maxrat)
c       print *,'R,BA,TA=',r_predef(iaa),ba_predef(iaa),ta_predef(iaa)
      end do
      nbranch=nbranch-1
      call trnsfi(nneig,nneig_sav,nslt)
      return
2000  format(' Number of atoms in the two branches exceed ',i5,
     -  ' - increase the parameter MAXBRANCH')
c2001  format(' Ring member:',i6,1x,a,1x,a,i5)
2002  format(' List of bonds in the ring:',/,
     -  (i3,' Bond',i6,1x,a,1x,a,i5,' - ',i6,1x,a,1x,a,i5))
2003  format(' Breaking bond ',i6,1x,a,1x,a,i5,' - ',i6,1x,a,1x,a,i5)
      end
