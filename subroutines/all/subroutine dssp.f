      subroutine dssp(c,n1,n,nslt,line,index,inamcol1,inamcol2,
     -  iresncol1,iresncol2,nneig,ineig,nbox,indices,ihbneig,ixc,ixo,
     -  ixn,ixa,dssplab,idistdssp,ch,cres,enghb,iparal,iantiparal,nss,
     -  itypss,ifss,ilss,ires0,nconfig,iwdssp,iwrose,iwhead,ifail,
     -  radtodeg,maxrepconf,maxng,nnlistlen,maxbox,listlen,maxss,
     -  maxrsd,maxrec)
      dimension c(3,n),index(n),nneig(n),ineig(maxng,nnlistlen),
     -  indices(maxbox,listlen),nbox(listlen),ihbneig(n),ixc(n),ixo(n),
     -  ixn(n),ixa(n),cres(3,n),iparal(n),iantiparal(n),
     -  idistdssp(9,maxrsd),ch(3,maxrsd),enghb(maxrsd),itypss(maxss),
     -  ifss(maxss),ilss(maxss)
      character*1 dssplab(maxrsd)
      dimension ccprev(3),rx(3),rn(3),rnprev(3),rn1(3),rn1prev(3)
      character*1 charl1(50),charl2(50)
      character*8 atnam
      character* 132 line(maxrec)
      character*1 typc
      character*21 ssname
      common /dsspnames/ lssname(9),ssname(9),typc(9)
      real*8 cosa
      data ing /0/,ils /0/
      lnam=inamcol2-inamcol1+1
c     print *,'DSSP n1,n,nslt,iwdssp=', n1,n,nslt,iwdssp
      ifail=0
      iok=1
      nres=0
      nresok=0
      iccprevf=0
      icccurrf=0
      icfound=0
      iofound=0
      infound=0
      ihfound=0
      iafound=0
      do ia=n1,nslt
        atnam(1:lnam)=line(index(ia))(inamcol1:inamcol2)
        if (lnam .eq. 4) atnam(5:8)='    '
        if (lnam .gt. 4) call leftadjustn(atnam,atnam,lnam)
        if (atnam(1:4) .eq. 'C   ' .or. atnam(1:4) .eq. ' C   ') then
          icfound=1
          ixc(nres+1)=ia
        else if (atnam(1:4) .eq. 'O   ' .or. atnam(1:4) .eq. ' O   ')
     -           then
          iofound=1
          ixo(nres+1)=ia
        else if (atnam(1:4) .eq. 'N   ' .or. atnam(1:4) .eq. ' N   ')
     -           then
          infound=1
          ixn(nres+1)=ia
        else if (atnam(1:4) .eq. 'CA  ' .or. atnam(1:4) .eq. ' CA  ')
     -           then
          iafound=1
          ixa(nres+1)=ia
        else if (atnam(1:4) .eq. 'H   ' .or. atnam(1:4) .eq. ' H   '
     -    .or. atnam(1:4) .eq. ' D  ' .or. atnam(1:4) .eq. 'D   ' .or.
     -    atnam(1:4) .eq. 'HN  ' .or. atnam(1:4) .eq. ' HN ') then
          ihfound=1
          call trnsfr(ch(1,nres+1),c(1,ia),3)
        end if
        if (ia .eq. nslt) then
          iresn=-1
        else
          read (line(index(ia+1))(iresncol1:iresncol2),*,ERR=999) iresn
          if (ia .eq. n1) ireso=iresn
        end if
        if (iresn .ne. ireso) then
c         New residue
          nres=nres+1
          iok=1
          if (ireso .ne. 0) then
c           print *,'IRESO,N=',ireso,iresn
c           print *,'ICFOUND,INFOUND,IOFOUND,IHFOUND=',
c    -          icfound,infound,iofound,ihfound
            if (icfound*infound*iofound .lt. 1) then
              if (nconfig .le. maxrepconf)
     -          print *,'Residue ',ireso,' is missing N, C or O'
              iok=0
            else if (ihfound .eq. 0) then
c             Generate H coordinates from C(prev), CA and N
              iok=0
              call zeroit(ch(1,nres),3)
              if (iccprevf*iafound .eq. 1) then
                dnh=0.0
                dnc=0.0
                do k=1,3
                  rx(k)=2.0*c(k,ixn(nres))-c(k,ixa(nres))-ccprev(k)
                  dnh=dnh+rx(k)**2
                  dnc=dnc+(c(k,ixn(nres))-ccprev(k))**2
                end do
                if (dnc .lt. 4.0) then
                  do k=1,3
                    ch(k,nres)=c(k,ixn(nres))+rx(k)/sqrt(dnh)
                  end do
                  iok=1
                else
c                 print *,'nconfig,maxrepconf=',nconfig,maxrepconf
                  if (nconfig .le. maxrepconf)
     -              print *,'Chain break at residue',nres,
     -              ' no H generated'
                end if
              else
                if (nconfig .le. maxrepconf)
     -            print *,'Could not generate H for residue',ireso
              end if
            end if
            if (icfound .eq. 1) then
              iccprevf=1
              if (nres .gt. 1 .and. ixn(nres) .gt. 0) then
                if (dist2(ccprev,c(1,ixn(nres))) .gt. 4.0) then
                  ihbneig(nres-1)=-1
                  if (nconfig .le. maxrepconf)
     -              print *,'Chain break between residues',nres-1,
     -                ' and',nres
                end if
              end if
              call trnsfr(ccprev,c(1,ixc(nres)),3)
            else
              iccprevf=0
              call zeroit(ccprev,3)
            end if
            if (iok .eq. 1) then
              nresok=nresok+1
              ihbneig(nres)=0
            else
              ihbneig(nres)=-1
            end if
          end if
          ireso=iresn
          icfound=0
          iofound=0
          infound=0
          ihfound=0
          iafound=0
        end if
c       write (77,*) 'ia,ireso,nres=',ia,ireso,nres
      end do
      if (nresok .eq. 0) then
        print *,'ERROR: No residues containing C, O and N were found'
        ifail=1
        return
      end if
c      do i=1,nres
c        write (77,1144) i,'CA:',line(index(ixa(i)))(1:80)
c        write (77,1144) i,'C :',line(index(ixc(i)))(1:80)
c        write (77,1144) i,'N :',line(index(ixn(i)))(1:80)
c        write (77,1144) i,'O :',line(index(ixo(i)))(1:80)
c        dd=dist2(ch(1,i),c(1,ixn(i)))
c        ddco=dist2(c(1,ixc(i)),c(1,ixo(i)))
c        ddcn=dist2(c(1,ixc(i)),c(1,ixn(i)))
c        write (77,*) 'dCO=',sqrt(ddco),' dCN=',sqrt(ddcn)
c        write (77,1133) i,sqrt(dd),(ch(k,i),k=1,3)
c1133    format(i5,' d(N-H)=',f10.5,' ch=',3f10.5)
c1144    format(i6,1x,a,1x,a)
c1145    format(i6,1x,a,1x,3f10.5)
c      end do
      threshold=-0.5/(0.42*0.20*332.0)
      do ir=1,nres
        enghb(ir)=threshold+1.0e-5
        dssplab(ir)=' '
        if (ixc(ir) .gt. 0) then
          call trnsfr(cres(1,ir),c(1,ixc(ir)),3)
        else
          call zeroit(cres(1,ir),3)
        end if
      end do
      rchb=9.2
      call nnlistsim(1,nres,cres,nneig,ineig,indices,nbox,
     -  rchb,ifail,maxng,nnlistlen,maxbox,listlen,0)
      if (ifail .gt. 0) then
        print *,'Linked-cell routine failed'
        do ir=1,nres
          if (ihbneig(ir) .ge. 0) then
            do jr=ir+3,nres
              if (ihbneig(jr) .ge. 0) then
                if (dist2(c(1,ixc(ir)),c(1,ixc(jr))) .lt. rchb**2) then
                  eij=1.0/sqrt(dist2(c(1,ixo(ir)),c(1,ixn(jr))))+
     -              1.0/sqrt(dist2(c(1,ixc(ir)),ch(1,jr)))-
     -              1.0/sqrt(dist2(c(1,ixo(ir)),ch(1,jr)))-
     -              1.0/sqrt(dist2(c(1,ixc(ir)),c(1,ixn(jr))))
                  eji=1.0/sqrt(dist2(c(1,ixo(jr)),c(1,ixn(ir))))+
     -              1.0/sqrt(dist2(c(1,ixc(jr)),ch(1,ir)))-
     -              1.0/sqrt(dist2(c(1,ixo(jr)),ch(1,ir)))-
     -              1.0/sqrt(dist2(c(1,ixc(jr)),c(1,ixn(ir))))
c                if (eij .lt. threshold .or. eji .lt. threshold)
c     -            write (77,1156) ir,jr,eij*(0.42*0.20*332.0),
c     -              eji*(0.42*0.20*332.0)
c1156            format(' ir,jr=',2i5,' eij,eji=',2f10.5)
                  if (eij .lt. enghb(ir)) then
                    if (enghb(ir) .eq. 0.0) write (6,1157) ir,jr,
     -                ihbneig(ir),enghb(ir),eij
                    ihbneig(ir)=jr
                    enghb(ir)=eij
                  end if
                  if (eji .lt. enghb(jr)) then
                    if (enghb(jr) .eq. 0.0) write (6,1157) jr,ir,
     -                ihbneig(jr),enghb(ir),eji
                    ihbneig(jr)=ir
                    enghb(jr)=eji
                  end if
                end if
              end if
            end do
          end if
        end do
      else
        do ir=1,nres
          if (ihbneig(ir) .ge. 0) then
            do jjr=1,nneig(ir)
              jr=ineig(jjr,ir)
              if (ihbneig(jr) .ge. 0 .and. jr .gt. ir+2) then
                if (dist2(c(1,ixc(ir)),c(1,ixc(jr))) .lt. rchb**2) then
                  eij=1.0/sqrt(dist2(c(1,ixo(ir)),c(1,ixn(jr))))+
     -              1.0/sqrt(dist2(c(1,ixc(ir)),ch(1,jr)))-
     -              1.0/sqrt(dist2(c(1,ixo(ir)),ch(1,jr)))-
     -              1.0/sqrt(dist2(c(1,ixc(ir)),c(1,ixn(jr))))
                  eji=1.0/sqrt(dist2(c(1,ixo(jr)),c(1,ixn(ir))))+
     -              1.0/sqrt(dist2(c(1,ixc(jr)),ch(1,ir)))-
     -              1.0/sqrt(dist2(c(1,ixo(jr)),ch(1,ir)))-
     -              1.0/sqrt(dist2(c(1,ixc(jr)),c(1,ixn(ir))))
c                  if (eij .lt. threshold .or. eji .lt. threshold)
c                  if (ir .lt. 999)
c     -              write (77,1155) ir,jr,eij,eji,threshold
c1155              format(' ir,jr=',2i5,' eij,eji=',2f10.5,' tr=',f10.5)
                  if (eij .lt. enghb(ir)) then
                    if (enghb(ir) .eq. 0.0) write (6,1157) ir,jr,
     -                ihbneig(ir),enghb(ir),eij
                    ihbneig(ir)=jr
                    enghb(ir)=eij
                  end if
                  if (eji .lt. enghb(jr)) then
                    if (enghb(jr) .eq. 0.0) write (6,1157) jr,ir,
     -                ihbneig(jr),enghb(ir),eji
                    ihbneig(jr)=ir
                    enghb(jr)=eji
                  end if
c                  if (ir .lt. 999) write (77,1154) ir,ihbneig(ir),
c     -              enghb(ir),jr,ihbneig(jr),enghb(jr)
c1154              format(' ir,ihbneig(ir),enghb(ir)=',2i5,f10.5,
c     -              ' jr,ihbneig(jr),enghb(jr)=',2i5,f10.5)
                end if
              end if
            end do
          end if
        end do
      end if
c     A hydrogen bond exist from C=O(i) to NH(ihbneig(i))
      nss=0
      call zeroiti(iparal,0,nres)
      call zeroiti(iantiparal,0,nres)
      ir=1
      do while (ir .lt. nres)
c       Look for the next SS element
        do while (ir .lt. nres .and. ihbneig(ir) .le. 0)
          ir=ir+1
        end do
        nhbinc=ihbneig(ir)-ir
        irf=ir
        nss0=nss
c       write (77,*) 'Start check ir=',ir,' nhbinc=',nhbinc
        if (nhbinc .gt. 2 .and. nhbinc .le. 5) then
c         Possible helix start
          ihfound=1
          nstep=1
          ir=ir+1
          do while (ihfound .eq. 1)
            nhbinc1=ihbneig(ir)-(ir)
            nhbinc2=ihbneig(ir+1)-(ir+1)
            nhbinc3=ihbneig(ir+2)-(ir+2)
            if (nhbinc1 .ne. nhbinc) then
              if (nhbinc2 .ne. nhbinc) then
                if (nhbinc1 .eq. nhbinc2) then
                  ihfound=0
                else if (nhbinc3 .ne. nhbinc) then
                  ihfound=0
                end if
              end if
            end if
c           write (77,6622) ir,nhbinc1,nhbinc2,nhbinc3,ihfound
c6622       format(' ir=',i4,' dihn1,2,3=',3i3,' ifhound=',i2)
            if (ihfound .eq. 1) then
              nstep=nstep+1
              irprev=ir
              ir=ir+1
            end if
          end do
          if (ir-irf .gt. 1) then
            nss=nss+1
            ifss(nss)=irf
            ilss(nss)=ir+nhbinc-1
            itypss(nss)=nhbinc+3
c           write (06,*) 'Helix nss,ifss,ilss,itypss=',nss,ifss(nss),
c    -        ilss(nss),itypss(nss)
          else
            ir=ir-1
          end if
        else if (nhbinc .eq. -3) then
c         Possible Lambda helix
          ihfound=1
          nstep=1
          ir=ir+1
          do while (ihfound .eq. 1)
            nhbinc1=ihbneig(ir)-(ir)
            if (nhbinc1 .ne. nhbinc) then
              ihfound=0
            end if
c           write (77,6622) ir,nhbinc1,nhbinc2,nhbinc3,ihfound
c6622       format(' ir=',i4,' dihn1,2,3=',3i3,' ifhound=',i2)
            if (ihfound .eq. 1) then
              nstep=nstep+1
              irprev=ir
              ir=ir+1
            end if
          end do
          if (ir-irf .gt. 1) then
            nss=nss+1
            ifss(nss)=irf
            ilss(nss)=ir+nhbinc-1
            itypss(nss)=9
c           write (77,*) 'L Helix nss,ifss,ilss,itypss=',nss,ifss(nss),
          end if
        else if (ihbneig(ir) .gt. 0) then
c         Possible sheet
          isfound=1
          isfirst=0
          islast=0
          npar=0
          napar=0
          ineigmax=0
          ineigmin=nres
          ifs=-1
          do while (isfound .eq. 1)
            jr=ihbneig(ir)
c           write (77,*) 'Sheet? ir,jr,ihbneig(jr)=',ir,jr,ihbneig(jr)
            if (jr .gt. 0) then
              if (ihbneig(jr) .eq. ir) then
c               HB(i,j) = HB(j,i) => Antiparallel
                iantiparal(ir)=1
                ing=jr
                ifs=ir-1
                ils=ir+1
c               write (77,*) 'A1 ir=',ir
              end if
            end if
            if (jr .gt. 2) then
              if (ihbneig(jr-2) .gt. 0) then
                if (ihbneig(jr-2) .eq. ir) then
c                 HB(j-1,i)=HB(i,j+1) => Parallel
                  iparal(ir)=1
                  ing=jr
                  ifs=ir
                  ils=ir+2
c                 write (77,*) 'P2 ir=',ir
                end if
              end if
            end if
            if (ir .gt. 1 .and. ir .lt. nres) then
              jrm=ihbneig(ir-1)
              if (jrm .gt. 0) then
                if (ihbneig(jrm) .eq. ir+1) then
c                 HB(i-1,j)=HB(j,i+1) => Parallel
                  iparal(ir)=1
                  ing=jrm
                  ifs=ir
                  ils=ir+2
c                 write (77,*) 'P1 ir=',ir
                end if
              end if
            else
              jrm=0
            end if
            if (jrm .gt. 2) then
              if (ihbneig(jrm-2) .gt. 0) then
                if (ihbneig(jrm-2) .eq. ir+1) then
c                 HB(i-1,j+1)=HB(j-1,i+1) => Antiparallel
                  iantiparal(ir)=1
                  ing=jrm
                  ifs=ir
                  ils=ir+2
c                 write (77,*) 'A2 ir=',ir
                end if
              end if
            end if
            if (iparal(ir)+iantiparal(ir) .gt. 0) then
              if (isfirst .eq. 0) isfirst=ifs
              if (ineigmin .gt. ing) ineigmin=ing
              if (ineigmax .lt. ing) ineigmax=ing
            else if (ir .eq. 1) then
              isfound=0
            else if (ir .gt. 1) then
              if (iparal(ir-1)+iantiparal(ir-1) .eq. 0) then
                isfound=0
                islast=ils
              end if
            end if
            npar=npar+iparal(ir)
            napar=napar+iantiparal(ir)
            ir=ir+1
c           write (77,*)'End do isfound,ir=',isfound,ir
          end do
c         write (77,*)'islast,isfirst,npar,napar,ineigmin,max=',
c    -       islast,isfirst,npar,napar,ineigmin,ineigmax
          if (islast-isfirst .gt. 2 .and. isfirst .gt. 0) then
            nss=nss+1
            ifss(nss)=isfirst
            ilss(nss)=islast
            if (npar*napar .gt. 0) then
              itypss(nss)=3
            else if (npar .gt. 0)  then
              itypss(nss)=1
              if (ineigmax-ineigmin .gt. (islast-isfirst)*2)
     -          itypss(nss)=4
            else if (napar .gt. 0)  then
              itypss(nss)=2
              if (ineigmax-ineigmin .gt. (islast-isfirst)*2)
     -          itypss(nss)=5
            end if
c           write (77,*) 'Sheet nss,ifss,ilss,itypss=',nss,ifss(nss),
c    -        ilss(nss),itypss(nss)
          end if
        end if
        if (nss .eq. maxss) then
          write (6,1002) maxss
          if (iwdssp .gt. 0) write (iwdssp,1002) maxss
          ir=nres
        else if (nss .eq. nss0) then
c         Neither helix nor sheet - just skip over
c         do while (ir .lt. nres .and. ihbneig(ir) .ne. 0)
            ir=ir+1
c         end do
        else
          if (nss .gt. 1) then
            if (ilss(nss-1) .ge. ifss(nss)) ifss(nss)=ilss(nss-1)+1
          else
            if (ifss(1) .lt. 1) ifss(1)=1
          end if
          ir=ir+1
        end if
      end do
      do iss=1,nss
c        write (77,1001) iss,ifss(iss),ilss(iss),itypss(iss)
c1001    format(' SS',i4,' Start at',i4,' End at',i4,' Type=',i2)
        do ir=ifss(iss),ilss(iss)
          dssplab(ir)=typc(itypss(iss))
          idistdssp(itypss(iss),ir)=idistdssp(itypss(iss),ir)+1
        end do
      end do
      if (iwdssp .gt. 0) then
        if (nconfig .le. 1) then
          if (nss .gt. 0) then
            write (iwdssp,2005) (i,ires0+ifss(i),ires0+ilss(i),
     -        ssname(itypss(i))(1:lssname(itypss(i))),i=1,nss)
            write (6,2005) (i,ires0+ifss(i),ires0+ilss(i),
     -        ssname(itypss(i))(1:lssname(itypss(i))),i=1,nss)
          else
            write (iwdssp,2006)
            write (6,2006)
          end if
        end if
        if (iwhead .eq. 1)
     -     write (iwdssp,2004) (typc(i),ssname(i)(1:lssname(i)),i=1,9)
        iresf=1
        do while (iresf .le. nres)
          iresl=min0(nres,iresf+49)
          do ic=1,50
            charl1(ic)=' '
            charl2(ic)=' '
          end do
          do ir=max0(3,iresf),min0(iresl,nres-2)
            if (ixa(ir-2)*ixa(ir)*ixa(ir+2) .ne. 0) then
              call angles(dist2(c(1,ixa(ir-2)),c(1,ixa(ir))),
     -          dist2(c(1,ixa(ir+2)),c(1,ixa(ir))),
     -          dist2(c(1,ixa(ir-2)),c(1,ixa(ir+2))),ca1,ca2,cbend)
              cosa=dble(cbend)
              bend=180.0-dacoscheck(cosa,ccc,1,6,'DSSP')*radtodeg
c             write (77,*) ir,' bend=',bend
              if (bend .gt. 70.0) charl1(ir-iresf+1)='S'
            else
              charl1(ir-iresf+1)='?'
            end if
          end do
          do ir=max0(2,iresf),min0(iresl,nres-2)
            if (ixa(ir-1)*ixa(ir)*ixa(ir+1)*ixa(ir+2) .ne. 0) then
              tors=dihangl(c,ixa(ir-1),ixa(ir),ixa(ir+1),ixa(ir+2),0,n)
c             write (77,*) ir,' tors=',tors*radtodeg
              if (tors .ge. 0.0) then
                charl2(ir-iresf+1)='+'
              else
                charl2(ir-iresf+1)='-'
              end if
            else
              charl2(ir-iresf+1)='?'
            end if
          end do
          write (iwdssp,2000) iresf,iresl,(mod(i,10),i=iresf,iresl)
          write (iwdssp,2001) (mod(ihbneig(i)/1000,10),i=iresf,iresl)
          write (iwdssp,2001) (mod(ihbneig(i)/100,10),i=iresf,iresl)
          write (iwdssp,2001) (mod(ihbneig(i)/10,10),i=iresf,iresl)
          write (iwdssp,2001) (mod(ihbneig(i),10),i=iresf,iresl)
          write (iwdssp,2002) (charl1(i-iresf+1),i=iresf,iresl)
          write (iwdssp,2002) (charl2(i-iresf+1),i=iresf,iresl)
          write (iwdssp,2002) (dssplab(i),i=iresf,iresl)
          iresf=iresl+1
        end do
      end if
      if (iwrose .gt. 0) then
c       Tentative alternative for turn detection (see Rose's 1977 paper)
        write (iwrose,*) 'Data for turn detection with GW Rose method'
        call normplane(c(1,ixa(1)),c(1,ixa(3)),c(1,ixa(5)),rnprev)
        call normplane(c(1,ixa(2)),c(1,ixa(3)),c(1,ixa(4)),rn1prev)
        do ir=4,nres-2
          call radcirc(c(1,ixa(ir-2)),c(1,ixa(ir)),c(1,ixa(ir+2)),r)
          call normplane(c(1,ixa(ir-2)),c(1,ixa(ir)),c(1,ixa(ir+2)),rn)
          rnn=scprod(rn,rnprev)
          call normplane(c(1,ixa(ir-1)),c(1,ixa(ir)),c(1,ixa(ir+1)),rn1)
          rnn1=scprod(rn1,rn1prev)
          call angles(dist2(c(1,ixa(ir-2)),c(1,ixa(ir))),
     -      dist2(c(1,ixa(ir+2)),c(1,ixa(ir))),
     -      dist2(c(1,ixa(ir-2)),c(1,ixa(ir+2))),ca1,ca2,cbend)
          cosa=dble(cbend)
          bend=180.0-dacoscheck(cosa,ccc,1,6,'DSSP')*radtodeg
          charl1(1)=' '
          if (bend .gt. 70.0) charl1(1)='S'
          write (iwrose,2003) ir,charl1(1),dssplab(ir),r,rnn,rnn1
          call trnsfr(rnprev,rn,3)
          call trnsfr(rn1prev,rn1,3)
        end do
      end if
      return
999   write (6,1000) ia,line
      if (iwdssp .gt. 0) write (iwdssp,1000) ia,line
      return
1000  format(' ERROR: invalid residue number for atom ',i6,':',/,a)
1002  format(' ERROR: maximum number of secondary structure elements (',
     -  i3,') has been reached',/,' - redimension the program ')
1157  format(' eHB update ',i4,' to ',i4,' eold=',f6.2,' old partner=',
     -  i4,' enew=',f6.2)
2000  format(/,i6,'-',i5,': ',50i1)
2001  format(14x,50i1)
2002  format(14x,50a1)
2003  format(i5,' DSSP labels:',a1,1x,a1,' rc=',f8.3,
     -  ' rn.rnprv(2)=',f8.4,' rn.rnprv(1)=',f8.4)
2004  format(/' line 1: residue number (mod 10)',/,
     -  ' lines 2-5: the digits of the residue number to which the ',/,
     -  12x,'residue number of line 1 is H-bonded (if any)',/,
     -  ' line 6: S for residues with bend angle ',
     -  '(CA(ir-2)-CA(ir)-CA(ir+2)) > 70 deg',/,
     -  ' line 7: + or -, the sign of the ',
     -  '(CA(ir-1)-CA(ir)-CA(ir+1)-CA(ir+2)) angle',/,
     -  ' line 8: secondary structure element type:',/,
     -  9(9x,a1,': ',a))
2005  format(/,(' SS#',i3,' Residue index range: [',i5,',',i5,'] Type:',
     -  a))
2006  format(' No secondary structure element was found')
      end
