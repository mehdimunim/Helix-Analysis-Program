      subroutine retroinverso(c,iat1,iat2,nslt,iatnum,line,index,chain,
     -  nsegm,iaincr,inamcol1,inamcol2,iresncol1,iresncol2,irescol1,
     -  irescol2,pi,ifail,nneig,ineig,mask,iout,maxng,maxrec)
      dimension c(3,nslt),iatnum(nslt),index(nslt),nneig(nslt),
     -  ineig(maxng,nslt),mask(nslt)
      dimension cn(3),ca(3),cb(3),cd(3),cg(3),cm_na(3),cm_bd(3),
     -  e(3),rh(3),cc(3),ch(3),chn(3),ch23(3),con(3),co1(3),cco(3),
     -  oco(3),ccc(3),ch1(3),ch2(3),ch3(3),cnn(3),cnh(3),co(3)
      character*1 chain,ans
      character*3 resnam3
      character*4 atnam4
      character*8 atnam,resnam,resnamn,hnname
      character* 132 line(maxrec)
c     print *,'RETROINVERSO nslt=',nslt,' iout=',iout
      call quiz(ans,itertyp,' ',' ',0,'terminus type',13,0,5,6,00)
      lnam=inamcol2-inamcol1+1
      ifail=0
      nres=0
      nresok=0
      icfound=0
      iofound=0
      infound=0
      ihfound=0
      iafound=0
      inf=0
      iaf=0
      ia2=iat2
      frocc=1.0
      bfac=0.0
      rco=1.4
      rcn=1.3
      rcc=1.6
      rnh=0.96
      rch=1.1
      resnam3='   '
      resnamn='        '
      if (chain .eq. ' ') chain='A'
      call indexit(mask,1,nslt,0)
      npro=0
      nwr=0
      ia=1
      icharmm=0
      iamber=0
      iambinc=1
      read (line(index(1))(iresncol1:iresncol2),*,ERR=999) iresninc
      if (itertyp .eq. 3) iresninc=iresninc-1
      iresninc=iresninc-1
c     print *,'IAT1,2=',iat1,iat2
      do ia=iat2,iat1,-1
        atnam(1:lnam)=line(index(ia))(inamcol1:inamcol2)
        call leftadjust4(atnam(1:4),atnam4)
c       print *,'IA=',ia,' ATNAM=',atnam(1:lnam)
        if (lnam .eq. 4) atnam(5:8)='    '
        if (lnam .gt. 4) call leftadjustn(atnam,atnam,lnam)
        if (atnam4 .eq. 'C   ' .or. atnam4 .eq. ' C   ') then
          icfound=ia
          call trnsfr(cc,c(1,ia),3)
        else if (atnam4 .eq. 'O   ' .or. atnam4 .eq. 'O1   ' .or.
     -           atnam4 .eq. 'OT1  ') then
          iofound=ia
          call trnsfr(co,c(1,ia),3)
        else if (atnam4 .eq. 'N   ') then
          infound=ia
          call trnsfr(cn,c(1,ia),3)
        else if (atnam4 .eq. 'CA  ') then
          iafound=ia
          call trnsfr(ca,c(1,ia),3)
        else if (atnam4 .eq. 'H   ' .or. atnam4 .eq. 'HN   ' .or.
     -           atnam4 .eq. 'H1  ') then
          ihfound=ia
          hnname=atnam
          call trnsfr(ch,c(1,ia),3)
        end if
        if (ia .eq. iat1) then
          iresn=-1
        else
          read (line(index(ia-1))(iresncol1:iresncol2),*,ERR=999) iresn
c         print *,'IRESN=',iresn
          if (ia .eq. iat2) ireso=iresn
        end if
        if (iresn .ne. ireso) then
c         New residue
          nres=nres+1
c         print *,'NRES=',nres
          ia1=ia
c         Now make the transformation
          resnam(1:irescol2-irescol1+1)=
     -      line(index(ia1))(irescol1:irescol2)
          ndel_extra=0
          if (iresn .eq. -1) then
c           Remove extra hydrogens from N (if any)
            nhf=0
            do inn=1,nneig(infound)
              in=ineig(inn,infound)
              atnam(1:lnam)=line(index(in))(inamcol1:inamcol2)
              call leftadjust4(atnam(1:4),atnam4)
c             print *,'IN=',in,' atno=',iatnum(in),' nhf=',nhf
              if (atnam4 .eq. 'H2   ' .or. atnam4 .eq. 'H3  ') then
                mask(in)=0
                ndel_extra=ndel_extra+1
              end if
            end do
c           print *,'H clean ndel=',ndel_extra
          end if
          if (ihfound .eq. 0) then
           if (resnam(1:3) .eq. 'PRO') then
             if (nneig(infound) .ne. 3) then
               write (6,1004) infound,nneig(infound),
     -           (ineig(i,infound),i=1,nneig(infound))
               stop
             else
              call addatom(2,cg,c(1,ineig(2,infound)),
     -          c(1,ineig(1,infound)),cn,ch,rnh,0.0,0.0,0.0,1,pi,1,
     -          ifail)
             end if
           else
              print *,'NOTE: Residue ',ireso,' is missing the amide H'
              nht=0
              do iaa=ia1,ia2
                atnam(1:lnam)=line(index(iaa))(inamcol1:inamcol2)
                call leftadjust4(atnam(1:4),atnam4)
                if (atnam4 .eq. 'H1  ' .or. atnam4 .eq. 'H2  ' .or.
     -              atnam4 .eq. 'H3  ') nht=nht+1
              end do
              if (nht .gt. 0) then
                write (6,1002)
                ifail=1
                return
              end if
c             Generate the amide H
              if (iresn .eq. -1) then
                do k=1,3
                  ccnk=(cn(k)+ca(k))/2.0
                  co1(k)=2.0*ccnk-cc(k)
                end do
              else
                iaa=ia
                notfound=1
                do while (iaa .gt. 1 .and. notfound .eq. 1)
                  atnam(1:lnam)=line(index(iaa))(inamcol1:inamcol2)
                  call leftadjust4(atnam(1:4),atnam4)
                  if (atnam4 .eq. 'C   ') then
                    notfound=0
                    call trnsfr(co1,c(1,iaa),3)
                  end if 
                  iaa=iaa-1
                end do
              end if
              call addatom(2,cg,co1,ca,cn,ch,rnh,0.0,0.0,0.0,1,pi,1,
     -          ifail)
            end if
          else
            mask(ihfound)=0
          end if
c         print *,'NRES=',nres,' IA1,2=',ia1,ia2
          if (nres .eq. 1) then
c           Remove extra oxygen from C (if any)
            do inn=1,nneig(icfound)
              in=ineig(inn,icfound)
              atnam(1:lnam)=line(index(in))(inamcol1:inamcol2)
              call leftadjust4(atnam(1:4),atnam4)
c             print *,'IN=',in,' name=',atnam(1:4)
              if (atnam4 .eq. 'O2  ' .or. atnam4 .eq. 'OT2 ') then
                mask(in)=0
                ndel_extra=ndel_extra+1
              end if
            end do
c           Remove ACE terminus atoms (if any)
            do i=ia1,ia2
              atnam(1:lnam)=line(index(i))(inamcol1:inamcol2)
              call leftadjust4(atnam(1:4),atnam4)
              if (atnam4 .eq. '1HH3' .or. atnam4 .eq. '1HH3' .or.
     -            atnam4 .eq. '3HH3' .or. atnam4 .eq. 'CH3 ' .or.
     -            atnam4 .eq. 'HT1 ' .or. atnam4 .eq. 'HT2 ' .or.
     -            atnam4 .eq. 'HT3 ' .or. atnam4 .eq. 'CAT ' .or.
     -            atnam4 .eq. 'NT  ' .or. atnam4 .eq. 'HNT ') then
                mask(i)=0
                ndel_extra=ndel_extra+1
              end if
              if (atnam4 .eq. 'CAT ') icharmm=1
              if (atnam4 .eq. 'CH3 ') iamber=1
            end do
            if (iamber  .eq. 1) then
              iambinc=0
              do i=ia1,ia2
                atnam(1:lnam)=line(index(i))(inamcol1:inamcol2)
                call leftadjust4(atnam(1:4),atnam4)
                if (atnam4 .eq. 'N   ' .or. atnam4 .eq. 'H   ') then
                  mask(i)=0
                  ndel_extra=ndel_extra+1
                end if
              end do
c             If Amber Nter was read, search for ca,ch,ch 
              if (iafound*icfound*iofound .eq. 0 .and.
     -            iamber .eq. 1) then
                i=iat2
                iresr=iresn
                do while (i .gt. iat1 .and. iresr .gt. iresn-1)
                  atnam(1:lnam)=line(index(i))(inamcol1:inamcol2)
                  call leftadjust4(atnam(1:4),atnam4)
                  if (atnam4 .eq. 'C   ') then
                    icfound=i
                    call trnsfr(cc,c(1,i),3)
                  else if (atnam4 .eq. 'O   ') then
                    iofound=i
                    call trnsfr(co,c(1,i),3)
                  else if (atnam4 .eq. 'CA  ') then
                    iafound=i
                    call trnsfr(ca,c(1,i),3)
                    resnam3=line(index(i))(irescol1:irescol1+2)
                  end if          
                  i=i-1
                  read (line(index(i))(iresncol1:iresncol2),*,
     -              ERR=999) iresr
                end do
c       print *,'CA=',ca
c       print *,'CO->H=',co
c       print *,'CC->N=',cc
c       print *,'RESNAM3=',resnam3
              end if
c       print *,'RESNAM3=',resnam3,' IAMBER=',iamber
            end if
            if (iamber .eq. 0) resnam3=resnam(1:3)
c           print *,'IA1,2=',ia1,ia2,' IAMBER=',iamber
c           print *,'O clean ndel=',ndel_extra
            if (itertyp .eq. 2) then
c             Zwitterion - generate H1,H2,H3
              ang=angleijk(c,nslt,iafound,icfound,iofound,6)*180.0/pi
              dang1=dihangl(c,iofound,icfound,iafound,infound,ihfound,
     -          maxrec)*180.0/pi
              dang2=dang1+120.0
              dang3=dang1-120.0
c             print *,'DANG123=',dang1,dang2,dang3
              call addatom(1,cg,cn,ca,cc,ch23,rnh,ang,dang2,0.0,1,pi,1,
     -          ifail)
              nwr=nwr+1
              write (iout,1000) iaincr+nwr,' H2 ',resnam3,chain,
     -          iresninc+nres,ch23,frocc,bfac
              call addatom(1,cg,cn,ca,cc,ch23,rnh,ang,dang3,0.0,1,pi,1,
     -          ifail)
              nwr=nwr+1
              write (iout,1000) iaincr+nwr,' H3 ',resnam3,chain,
     -          iresninc+nres,ch23,frocc,bfac
              ndel_extra=ndel_extra-2
              iresninc=iresninc-iamber
            else if (itertyp .eq. 3) then
c             Blocked - generate ACE as separate residue (Amber)
              call addatom(2,cg,ca,co,cc,cco,rcn,0.0,0.0,57.0,1,pi,1,
     -          ifail)
              call addatom(1,cg,ca,cc,cco,oco,rco,120.0,90.0,0.0,1,
     -          pi,1,ifail)
              call addatom(2,cg,oco,cc,cco,ccc,rcc,0.0,0.0,0.0,1,
     -          pi,1,ifail)
              call addatom(1,cg,cn,cco,ccc,ch1,rch,104.0,0.0,0.0,1,
     -          pi,1,ifail)
              call addatom(1,cg,cn,cco,ccc,ch2,rch,104.0,120.0,0.0,1,
     -          pi,1,ifail)
              call addatom(1,cg,cn,cco,ccc,ch3,rch,104.0,-120.0,0.0,1,
     -          pi,1,ifail)
              nwr=nwr+1
              write (iout,1000) iaincr+nwr,'1HH3','ACE',chain,
     -          iresninc+nres-iambinc,ch1,frocc,bfac
              nwr=nwr+1
              write (iout,1000) iaincr+nwr,'2HH3','ACE',chain,
     -          iresninc+nres-iambinc,ch2,frocc,bfac
              nwr=nwr+1
              write (iout,1000) iaincr+nwr,'3HH3','ACE',chain,
     -          iresninc+nres-iambinc,ch3,frocc,bfac
              nwr=nwr+1
              write (iout,1000) iaincr+nwr,' CH3','ACE',chain,
     -          iresninc+nres-iambinc,ccc,frocc,bfac
              nwr=nwr+1
              write (iout,1000) iaincr+nwr,' C  ','ACE',chain,
     -          iresninc+nres-iambinc,cco,frocc,bfac
              nwr=nwr+1
              write (iout,1000) iaincr+nwr,' O  ','ACE',chain,
     -          iresninc+nres-iambinc,oco,frocc,bfac
            else if (itertyp .eq. 4) then
c             Blocked - generate ACE as same residue
c       print *,'CA=',ca
c       print *,'CC=',cc
c       print *,'CO=',co
c       print *,'RESNAM3=',resnam3
              call addatom(2,cg,ca,co,cc,cco,rcn,0.0,0.0,57.0,1,pi,1,
     -          ifail)
              call addatom(1,cg,ca,cc,cco,oco,rco,120.0,90.0,0.0,1,
     -          pi,1,ifail)
              call addatom(2,cg,oco,cc,cco,ccc,rcc,0.0,0.0,0.0,1,
     -          pi,1,ifail)
              call addatom(1,cg,cn,cco,ccc,ch1,rch,104.0,0.0,0.0,1,
     -          pi,1,ifail)
              call addatom(1,cg,cn,cco,ccc,ch2,rch,104.0,120.0,0.0,1,
     -          pi,1,ifail)
              call addatom(1,cg,cn,cco,ccc,ch3,rch,104.0,-120.0,0.0,1,
     -          pi,1,ifail)
              nwr=nwr+1
              write (iout,1000) iaincr+nwr,' HY1',resnam3,chain,
     -          iresninc+nres,ch1,frocc,bfac
              nwr=nwr+1
              write (iout,1000) iaincr+nwr,' HY2',resnam3,chain,
     -          iresninc+nres,ch2,frocc,bfac
              nwr=nwr+1
              write (iout,1000) iaincr+nwr,' HY3',resnam3,chain,
     -          iresninc+nres,ch3,frocc,bfac
              nwr=nwr+1
              write (iout,1000) iaincr+nwr,' CAY',resnam3,chain,
     -          iresninc+nres,ccc,frocc,bfac
              nwr=nwr+1
              write (iout,1000) iaincr+nwr,' CY ',resnam3,chain,
     -          iresninc+nres,cco,frocc,bfac
              nwr=nwr+1
              write (iout,1000) iaincr+nwr,' OY ',resnam3,chain,
     -          iresninc+nres,oco,frocc,bfac
              if (iamber .eq. 1) iresninc=iresninc-1
            end if
          end if
          if (iresn .eq. -1) then
c           Remove N terminal H and/or NME
            do i=ia1,ia2
              atnam(1:lnam)=line(index(i))(inamcol1:inamcol2)
              call leftadjust4(atnam(1:4),atnam4)
              if (atnam4 .eq. '1HH3' .or. atnam4 .eq. '1HH3' .or.
     -            atnam4 .eq. '3HH3' .or. atnam4 .eq. 'CH3 ' .or.
     -            atnam4 .eq. 'HY1 ' .or. atnam4 .eq. 'HY2 ' .or.
     -            atnam4 .eq. 'HY3 ' .or. atnam4 .eq. 'CAY ' .or.
     -            atnam4 .eq. 'CY  ' .or. atnam4 .eq. 'OY  ') then
                mask(i)=0
                ndel_extra=ndel_extra+1
              end if
              if (atnam4 .eq. 'CAY ') icharmm=1
              if (atnam4 .eq. 'CH3 ') iamber=1
            end do
          end if
c         write (6,8932) ia1,ia2,nres,iresn,iamber
c8932     format(' IA1,2=',2i5,' NRES,IRESN=',2i5,' IAMBER=',i2)
          if (iamber .eq. 1) resnamn(1:irescol2-irescol1+1)=
     -      line(index(ia1-1))(irescol1:irescol2)
          if (iamber .eq. 0 .or. (nres .gt. 1 .and. iresn .ne. -1)) then
c           write (6,9861) nres,ia1,ia2,resnam(1:3),ireso,iresn,
c    -        icfound,infound,iofound,ihfound,iafound
c9861       format(' nres=',i4,' ia1,2=',2i6,' resnam=',a,' ireso,n=',
c    -        2i4,/,' icfound,infound,iofound,ihfound,iafound=',5i6)
            if (icfound*infound*iofound*iafound .eq. 0) then
               print *,'ERROR: Residue ',ireso,
     -         ' is missing N, C, CA or O'
               ifail=1
               return
            end if
            hnname=' N      '
            mask(icfound)=0
            mask(infound)=0
            mask(iafound)=0
            mask(iofound)=0
            if (resnam(1:3) .eq. 'PRO') then
c             Special procedure for proline
              npro=npro+1
              ibfound=0
              idfound=0
              igfound=0
              ihfound=0
c             print *,'PRO ia1,ia2=',ia1,ia2
              do i=ia1,ia2
                atnam(1:lnam)=line(index(i))(inamcol1:inamcol2)
                call leftadjust4(atnam(1:4),atnam4)
                if (atnam4 .eq. 'CB  ') then
                  ibfound=i
                else if (atnam4 .eq. 'CD  ') then
                  idfound=i
                else if (atnam4 .eq. 'CG  ') then
                  igfound=i
                else if (atnam4 .eq. 'HA  ') then
                  ihfound=i
                  call trnsfr(ch,c(1,ihfound),3)
                end if
              end do
              if (ibfound*idfound*igfound*ihfound .eq. 0) then
                print *,'Proline residue ',ireso,
     -            ' has no CB, CD, CG or HA'
                ifail=1
                return
              end if
              nwr=nwr+1
              write (iout,1000) iaincr+nwr,' N  ','PRO',chain,
     -          iresninc+nres,(c(k,icfound),k=1,3),frocc,bfac
              call trnsfr(cn,c(1,icfound),3)
              nwr=nwr+1
              write (iout,1000) iaincr+nwr,' CA ','PRO',chain,
     -          iresninc+nres,(c(k,iafound),k=1,3),frocc,bfac
              call trnsfr(ca,c(1,iafound),3)
              nwr=nwr+1
              write (iout,1000) iaincr+nwr,' CB ','PRO',chain,
     -          iresninc+nres,(c(k,ibfound),k=1,3),frocc,bfac
              call trnsfr(cb,c(1,ibfound),3)
              nwr=nwr+1
              write (iout,1000) iaincr+nwr,' CD ','PRO',chain,
     -          iresninc+nres,(c(k,iofound),k=1,3),frocc,bfac
              call trnsfr(cd,c(1,iofound),3)
c              write (6,9866) icfound,iafound,ibfound,iofound,infound
c9866          format(' icfound,iafound,ibfound,iofound,infound=',5i6)
c             write (6,*) 'cn=',cn
c             write (6,*) 'ca=',ca
c             write (6,*) 'cb=',cb
c             write (6,*) 'cd=',cd
              phirad=dihangl(c,iofound,icfound,iafound,ibfound,0,maxrec)
     -          *180.0/pi
c             print *,'phirad=',phirad,' ',phirad
c             Generate new CG
              if (abs(phirad) .gt. 45.0) then
                print *,'NOTE: Proline ring # ',iresninc+nres,
     -            ' will be wide open'
                call addatom(1,cd,cn,ca,cb,cg,1.0,105.0,phirad,0.0,0,pi,
     -            1,ifail)
              else
                do k=1,3
                  cm_na(k)=(c(k,infound)+ca(k))/2.0
                end do
c               write (6,*) 'cm_na=',cm_na
                dcg=sqrt(dist2(c(1,igfound),cm_na))
                do k=1,3
                  cm_na(k)=(cn(k)+ca(k))/2.0
                  cm_bd(k)=(cb(k)+cd(k))/2.0
                end do
c               write (6,*) 'cm_na=',cm_na
c               write (6,*) 'cm_bd=',cm_bd
                call arrdiff(cm_bd,cm_na,e,3)
c               write (6,*) 'e=',e
                enorm=sqrt(scprod(e,e))
c               write (6,*) 'dcg=',dcg,' enorm=',enorm
                do k=1,3
                  cg(k)=cm_na(k)+dcg*e(k)/enorm
                end do
              end if
c             write (6,*) 'cg=',cg
              nwr=nwr+1
              write (iout,1000) iaincr+nwr,' CG ','PRO',chain,
     -          iresninc+nres,cg,frocc,bfac
c             Generate new hydrogens
              call addatom(2,e,ca,cg,cb,rh,1.0,0.0,0.0,57.,0,pi,1,ifail)
              nwr=nwr+1
              write (iout,1000) iaincr+nwr,' HB2','PRO',chain,
     -          iresninc+nres,rh,frocc,bfac
              call addatom(2,e,ca,cg,cb,rh,1.0,0.0,0.,-57.,0,pi,1,ifail)
              nwr=nwr+1
              write (iout,1000) iaincr+nwr,' HB3','PRO',chain,
     -          iresninc+nres,rh,frocc,bfac
              call addatom(2,e,cb,cd,cg,rh,1.0,0.0,0.0,57.,0,pi,1,ifail)
              nwr=nwr+1
              write (iout,1000) iaincr+nwr,' HG2','PRO',chain,
     -          iresninc+nres,rh,frocc,bfac
              call addatom(2,e,cb,cd,cg,rh,1.0,0.0,0.,-57.,0,pi,1,ifail)
              nwr=nwr+1
              write (iout,1000) iaincr+nwr,' HG3','PRO',chain,
     -          iresninc+nres,rh,frocc,bfac
              call addatom(2,e,cn,cg,cd,rh,1.0,0.0,0.0,57.,0,pi,1,ifail)
              nwr=nwr+1
              write (iout,1000) iaincr+nwr,' HD2','PRO',chain,
     -          iresninc+nres,rh,frocc,bfac
              call addatom(2,e,cn,cg,cd,rh,1.0,0.0,0.,-57.,0,pi,1,ifail)
              nwr=nwr+1
              write (iout,1000) iaincr+nwr,' HD3','PRO',chain,
     -          iresninc+nres,rh,frocc,bfac
              nwr=nwr+1
              write (iout,1000) iaincr+nwr,' HA ','PRO',chain,
     -          iresninc+nres,(c(k,ihfound),k=1,3),frocc,bfac
              nwr=nwr+1
              write (iout,1000) iaincr+nwr,' C  ','PRO',chain,
     -          iresninc+nres,(c(k,infound),k=1,3),frocc,bfac
              nwr=nwr+1
              write (iout,1000) iaincr+nwr,' O  ','PRO',chain,
     -          iresninc+nres,(c(k,idfound),k=1,3),frocc,bfac
            else
c             Write new backbone
              nsc=ia2-ia1-4-ndel_extra
              rco=sqrt(dist2(c(1,icfound),c(1,iofound)))
              rnh=sqrt(dist2(c(1,infound),ch))
              nwr=nwr+1
              write (iout,1000) iaincr+nwr,' N  ',resnam(1:3),chain,
     -          iresninc+nres,(c(k,icfound),k=1,3),frocc,bfac
              call rescale_bl(c(1,icfound),c(1,iofound),chn,rco,rnh)
              nwr=nwr+1
              if (nres .eq. 1 .and. itertyp .eq. 2) then
                write (iout,1000) iaincr+nwr,' H1 ',resnam(1:3),chain,
     -            iresninc+nres,chn,frocc,bfac
              else
                write (iout,1000) iaincr+nwr,' H  ',resnam(1:3),chain,
     -            iresninc+nres,chn,frocc,bfac
              end if
              nwr=nwr+1
              write (iout,1000) iaincr+nwr,' CA ',resnam(1:3),chain,
     -          iresninc+nres,(c(k,iafound),k=1,3),frocc,bfac
              do i=ia1,ia2
                if (mask(i) .gt. 0) then
                  nwr=nwr+1
                  write (iout,1000) iaincr+nwr,
     -              line(index(i))(inamcol1:inamcol2),resnam(1:3),chain,
     -              iresninc+nres,(c(k,i),k=1,3),frocc,bfac
                end if
              end do
              write (iout,1000) iaincr+nwr,' C  ',resnam(1:3),chain,
     -          iresninc+nres,(c(k,infound),k=1,3),frocc,bfac
              call rescale_bl(c(1,infound),ch,con,rnh,rco)
              nwr=nwr+1
              if ((iresn .eq. -1 .or. resnamn(1:3) .eq. 'ACE') .and.
     -             itertyp .eq. 2) then
                write (iout,1000) iaincr+nwr,' OT1',resnam(1:3),chain,
     -            iresninc+nres,con,frocc,bfac
              else
                write (iout,1000) iaincr+nwr,' O  ',resnam(1:3),chain,
     -          iresninc+nres,con,frocc,bfac
              end if
            end if
            inf=infound
            iaf=iafound
            icfound=0
            iofound=0
            infound=0
            ihfound=0
            iafound=0
            resnam3=resnam(1:3)
          end if
          ireso=iresn
c         print *,'IA1,IA2,IA,NRES=',ia1,ia2,ia,nres
          ia2=ia1-1
        end if
      end do
      if (itertyp .eq. 2) then
c       Zwitterion - generate O2
        call addatom(2,cg,con,c(1,iaf),c(1,inf),e,rco,0.0,0.0,
     -    0.0,0,pi,1,ifail)
        nwr=nwr+1
        write (iout,1000) iaincr+nwr,' OT2',resnam3,chain,
     -    iresninc+nres-iamber,e,frocc,bfac
      else if (itertyp .eq. 3) then
c       Blocked - generate NME as separate residue
c       print *,'CA=',ca
c       print *,'CN=',cn
c       print *,'CH=',ch
        call addatom(2,cg,ca,con,cn,cnn,rcn,0.0,0.0,0.0,1,pi,1,ifail)
        call addatom(1,cg,ca,cn,cnn,cnh,rnh,104.0,-120.0,0.0,1,
     -          pi,1,ifail)
        call addatom(2,cg,cn,cnh,cnn,ccc,rcn,0.0,0.0,0.0,1,pi,1,ifail)
        call addatom(1,cg,cn,cnn,ccc,ch1,rch,104.180,0.0,0.0,1,pi,1,
     -    ifail)
        call addatom(1,cg,cn,cnn,ccc,ch2,rch,104.0,120.0,0.0,1,pi,1,
     -    ifail)
        call addatom(1,cg,cn,cnn,ccc,ch3,rch,104.0,-120.0,0.0,1,pi,1,
     -    ifail)
        nwr=nwr+1
        write (iout,1000) iaincr+nwr,' N  ','NME',chain,
     -    iresninc+nres+iambinc,cnn,frocc,bfac
        nwr=nwr+1
        write (iout,1000) iaincr+nwr,' H  ','NME',chain,
     -    iresninc+nres+iambinc,cnh,frocc,bfac
        nwr=nwr+1
        write (iout,1000) iaincr+nwr,' CH3','NME',chain,
     -    iresninc+nres+iambinc,ccc,frocc,bfac
        nwr=nwr+1
        write (iout,1000) iaincr+nwr,'1HH3','NME',chain,
     -    iresninc+nres+iambinc,ch1,frocc,bfac
        nwr=nwr+1
        write (iout,1000) iaincr+nwr,'2HH3','NME',chain,
     -    iresninc+nres+iambinc,ch2,frocc,bfac
        nwr=nwr+1
        write (iout,1000) iaincr+nwr,'3HH3','NME',chain,
     -    iresninc+nres+iambinc,ch3,frocc,bfac
      else if (itertyp .eq. 4) then
c       Blocked - generate NME as same residue
        if (iamber .eq. 0) resnam3=resnam
        call addatom(2,cg,ca,con,cn,cnn,rcn,0.0,0.0,0.0,1,pi,1,ifail)
        call addatom(1,cg,ca,cn,cnn,cnh,rnh,104.0,-120.0,0.0,1,
     -          pi,1,ifail)
        call addatom(2,cg,cn,cnh,cnn,ccc,rcn,0.0,0.0,0.0,1,pi,1,ifail)
        call addatom(1,cg,cn,cnn,ccc,ch1,rch,104.0,0.0,0.0,1,pi,1,ifail)
        call addatom(1,cg,cn,cnn,ccc,ch2,rch,104.0,120.0,0.0,1,pi,1,
     -    ifail)
        call addatom(1,cg,cn,cnn,ccc,ch3,rch,104.0,-120.0,0.0,1,pi,1,
     -    ifail)
        nwr=nwr+1
        write (iout,1000) iaincr+nwr,' NT ',resnam3,chain,
     -    iresninc+nres-iamber,cnn,frocc,bfac
        nwr=nwr+1
        write (iout,1000) iaincr+nwr,' HNT',resnam3,chain,
     -    iresninc+nres-iamber,cnh,frocc,bfac
        nwr=nwr+1
        write (iout,1000) iaincr+nwr,' CAT',resnam3,chain,
     -    iresninc+nres-iamber,ccc,frocc,bfac
        nwr=nwr+1
        write (iout,1000) iaincr+nwr,' HT1',resnam3,chain,
     -    iresninc+nres-iamber,ch1,frocc,bfac
        nwr=nwr+1
        write (iout,1000) iaincr+nwr,' HT2',resnam3,chain,
     -    iresninc+nres-iamber,ch2,frocc,bfac
        nwr=nwr+1
        write (iout,1000) iaincr+nwr,' HT3',resnam3,chain,
     -    iresninc+nres-iamber,ch3,frocc,bfac
      end if
      if (npro .gt. 0) write (6,1001) npro
      return
999   print *,'Illegal record in line'
      print *,line(index(ia))
1000  format('ATOM  ',i5,1x,a4,1x,a3,1x,a1,i4,1x,3x,3f8.3,2f6.2)
1001  format(' NOTE:',i3,' prolines were found - run a minimzation to ',
     -  'correct the distorsions')
1002  format(' It looks like the input is not a single peptide',/,
     -  ' Run each peptide (chain) sperately')
1003  format(' Use the Modif<Y> option to add amide hydrogens')
1004  format(' PRO N (ia=',i7,') has no H but the # of bonds is not 3 ',
     -  'but',i3,/,' in:',10i7)
      end
