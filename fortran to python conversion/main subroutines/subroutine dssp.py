      def dssp(c,n1,n,nslt,line,index,inamcol1,inamcol2,
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
      #real*8 cosa
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
      for ia in range(n1, nslt):
        atnam(1:lnam)=line(index(ia))(inamcol1:inamcol2)
        if (lnam  ==  4) atnam(5:8)='    '
        if (lnam  >  4) call leftadjustn(atnam,atnam,lnam)
        if (atnam(1:4)  ==  'C   '  or  atnam(1:4)  ==  ' C   ') :
          icfound=1
          ixc(nres+1)=ia
        elif (atnam(1:4)  ==  'O   '  or  atnam(1:4)  ==  ' O   ')
     -           :
          iofound=1
          ixo(nres+1)=ia
        elif (atnam(1:4)  ==  'N   '  or  atnam(1:4)  ==  ' N   ')
     -           :
          infound=1
          ixn(nres+1)=ia
        elif (atnam(1:4)  ==  'CA  '  or  atnam(1:4)  ==  ' CA  ')
     -           :
          iafound=1
          ixa(nres+1)=ia
        elif (atnam(1:4)  ==  'H   '  or  atnam(1:4)  ==  ' H   '
     -     or  atnam(1:4)  ==  ' D  '  or  atnam(1:4)  ==  'D   '  or 
     -    atnam(1:4)  ==  'HN  '  or  atnam(1:4)  ==  ' HN ') :
          ihfound=1
          call trnsfr(ch(1,nres+1),c(1,ia),3)
        ## end if
        if (ia  ==  nslt) :
          iresn=-1
        else:
          read (line(index(ia+1))(iresncol1:iresncol2),*,ERR=999) iresn
          if (ia  ==  n1) ireso=iresn
        ## end if
        if (iresn  !=  ireso) :
c         New residue
          nres=nres+1
          iok=1
          if (ireso  !=  0) :
c           print *,'IRESO,N=',ireso,iresn
c           print *,'ICFOUND,INFOUND,IOFOUND,IHFOUND=',
c    -          icfound,infound,iofound,ihfound
            if (icfound*infound*iofound  <  1) :
              if (nconfig  <=  maxrepconf)
     -          print *,'Residue ',ireso,' is missing N, C or O'
              iok=0
            elif (ihfound  ==  0) :
c             Generate H coordinates from C(prev), CA and N
              iok=0
              call zeroit(ch(1,nres),3)
              if (iccprevf*iafound  ==  1) :
                dnh=0.0
                dnc=0.0
                for k in range(0, 3):
                  rx(k)=2.0*c(k,ixn(nres))-c(k,ixa(nres))-ccprev(k)
                  dnh=dnh+rx(k)**2
                  dnc=dnc+(c(k,ixn(nres))-ccprev(k))**2
                
                if (dnc  <  4.0) :
                  for k in range(0, 3):
                    ch(k,nres)=c(k,ixn(nres))+rx(k)/sqrt(dnh)
                  
                  iok=1
                else:
c                 print *,'nconfig,maxrepconf=',nconfig,maxrepconf
                  if (nconfig  <=  maxrepconf)
     -              print *,'Chain break at residue',nres,
     -              ' no H generated'
                ## end if
              else:
                if (nconfig  <=  maxrepconf)
     -            print *,'Could not generate H for residue',ireso
              ## end if
            ## end if
            if (icfound  ==  1) :
              iccprevf=1
              if (nres  >  1  and  ixn(nres)  >  0) :
                if (dist2(ccprev,c(1,ixn(nres)))  >  4.0) :
                  ihbneig(nres-1)=-1
                  if (nconfig  <=  maxrepconf)
     -              print *,'Chain break between residues',nres-1,
     -                ' and',nres
                ## end if
              ## end if
              call trnsfr(ccprev,c(1,ixc(nres)),3)
            else:
              iccprevf=0
              call zeroit(ccprev,3)
            ## end if
            if (iok  ==  1) :
              nresok=nresok+1
              ihbneig(nres)=0
            else:
              ihbneig(nres)=-1
            ## end if
          ## end if
          ireso=iresn
          icfound=0
          iofound=0
          infound=0
          ihfound=0
          iafound=0
        ## end if
c       write (77,*) 'ia,ireso,nres=',ia,ireso,nres
      
      if (nresok  ==  0) :
        print *,'ERROR: No residues containing C, O and N were found'
        ifail=1
        return
      ## end if
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
c      
      threshold=-0.5/(0.42*0.20*332.0)
      for ir in range(0, nres):
        enghb(ir)=threshold+1.0e-5
        dssplab(ir)=' '
        if (ixc(ir)  >  0) :
          call trnsfr(cres(1,ir),c(1,ixc(ir)),3)
        else:
          call zeroit(cres(1,ir),3)
        ## end if
      
      rchb=9.2
      call nnlistsim(1,nres,cres,nneig,ineig,indices,nbox,
     -  rchb,ifail,maxng,nnlistlen,maxbox,listlen,0)
      if (ifail  >  0) :
        print *,'Linked-cell routine failed'
        for ir in range(0, nres):
          if (ihbneig(ir)  >=  0) :
            for jr in range(ir+3, nres):
              if (ihbneig(jr)  >=  0) :
                if (dist2(c(1,ixc(ir)),c(1,ixc(jr)))  <  rchb**2) :
                  eij=1.0/sqrt(dist2(c(1,ixo(ir)),c(1,ixn(jr))))+
     -              1.0/sqrt(dist2(c(1,ixc(ir)),ch(1,jr)))-
     -              1.0/sqrt(dist2(c(1,ixo(ir)),ch(1,jr)))-
     -              1.0/sqrt(dist2(c(1,ixc(ir)),c(1,ixn(jr))))
                  eji=1.0/sqrt(dist2(c(1,ixo(jr)),c(1,ixn(ir))))+
     -              1.0/sqrt(dist2(c(1,ixc(jr)),ch(1,ir)))-
     -              1.0/sqrt(dist2(c(1,ixo(jr)),ch(1,ir)))-
     -              1.0/sqrt(dist2(c(1,ixc(jr)),c(1,ixn(ir))))
c                if (eij  <  threshold  or  eji  <  threshold)
c     -            write (77,1156) ir,jr,eij*(0.42*0.20*332.0),
c     -              eji*(0.42*0.20*332.0)
c1156            format(' ir,jr=',2i5,' eij,eji=',2f10.5)
                  if (eij  <  enghb(ir)) :
                    if (enghb(ir)  ==  0.0) write (6,1157) ir,jr,
     -                ihbneig(ir),enghb(ir),eij
                    ihbneig(ir)=jr
                    enghb(ir)=eij
                  ## end if
                  if (eji  <  enghb(jr)) :
                    if (enghb(jr)  ==  0.0) write (6,1157) jr,ir,
     -                ihbneig(jr),enghb(ir),eji
                    ihbneig(jr)=ir
                    enghb(jr)=eji
                  ## end if
                ## end if
              ## end if
            
          ## end if
        
      else:
        for ir in range(0, nres):
          if (ihbneig(ir)  >=  0) :
            for jjr in range(0, nneig(ir)):
              jr=ineig(jjr,ir)
              if (ihbneig(jr)  >=  0  and  jr  >  ir+2) :
                if (dist2(c(1,ixc(ir)),c(1,ixc(jr)))  <  rchb**2) :
                  eij=1.0/sqrt(dist2(c(1,ixo(ir)),c(1,ixn(jr))))+
     -              1.0/sqrt(dist2(c(1,ixc(ir)),ch(1,jr)))-
     -              1.0/sqrt(dist2(c(1,ixo(ir)),ch(1,jr)))-
     -              1.0/sqrt(dist2(c(1,ixc(ir)),c(1,ixn(jr))))
                  eji=1.0/sqrt(dist2(c(1,ixo(jr)),c(1,ixn(ir))))+
     -              1.0/sqrt(dist2(c(1,ixc(jr)),ch(1,ir)))-
     -              1.0/sqrt(dist2(c(1,ixo(jr)),ch(1,ir)))-
     -              1.0/sqrt(dist2(c(1,ixc(jr)),c(1,ixn(ir))))
c                  if (eij  <  threshold  or  eji  <  threshold)
c                  if (ir  <  999)
c     -              write (77,1155) ir,jr,eij,eji,threshold
c1155              format(' ir,jr=',2i5,' eij,eji=',2f10.5,' tr=',f10.5)
                  if (eij  <  enghb(ir)) :
                    if (enghb(ir)  ==  0.0) write (6,1157) ir,jr,
     -                ihbneig(ir),enghb(ir),eij
                    ihbneig(ir)=jr
                    enghb(ir)=eij
                  ## end if
                  if (eji  <  enghb(jr)) :
                    if (enghb(jr)  ==  0.0) write (6,1157) jr,ir,
     -                ihbneig(jr),enghb(ir),eji
                    ihbneig(jr)=ir
                    enghb(jr)=eji
                  ## end if
c                  if (ir  <  999) write (77,1154) ir,ihbneig(ir),
c     -              enghb(ir),jr,ihbneig(jr),enghb(jr)
c1154              format(' ir,ihbneig(ir),enghb(ir)=',2i5,f10.5,
c     -              ' jr,ihbneig(jr),enghb(jr)=',2i5,f10.5)
                ## end if
              ## end if
            
          ## end if
        
      ## end if
c     A hydrogen bond exist from C=O(i) to NH(ihbneig(i))
      nss=0
      call zeroiti(iparal,0,nres)
      call zeroiti(iantiparal,0,nres)
      ir=1
      while (ir  <  nres)
c       Look for the next SS element
        while (ir  <  nres  and  ihbneig(ir)  <=  0)
          ir=ir+1
        
        nhbinc=ihbneig(ir)-ir
        irf=ir
        nss0=nss
c       write (77,*) 'Start check ir=',ir,' nhbinc=',nhbinc
        if (nhbinc  >  2  and  nhbinc  <=  5) :
c         Possible helix start
          ihfound=1
          nstep=1
          ir=ir+1
          while (ihfound  ==  1)
            nhbinc1=ihbneig(ir)-(ir)
            nhbinc2=ihbneig(ir+1)-(ir+1)
            nhbinc3=ihbneig(ir+2)-(ir+2)
            if (nhbinc1  !=  nhbinc) :
              if (nhbinc2  !=  nhbinc) :
                if (nhbinc1  ==  nhbinc2) :
                  ihfound=0
                elif (nhbinc3  !=  nhbinc) :
                  ihfound=0
                ## end if
              ## end if
            ## end if
c           write (77,6622) ir,nhbinc1,nhbinc2,nhbinc3,ihfound
c6622       format(' ir=',i4,' dihn1,2,3=',3i3,' ifhound=',i2)
            if (ihfound  ==  1) :
              nstep=nstep+1
              irprev=ir
              ir=ir+1
            ## end if
          
          if (ir-irf  >  1) :
            nss=nss+1
            ifss(nss)=irf
            ilss(nss)=ir+nhbinc-1
            itypss(nss)=nhbinc+3
c           write (06,*) 'Helix nss,ifss,ilss,itypss=',nss,ifss(nss),
c    -        ilss(nss),itypss(nss)
          else:
            ir=ir-1
          ## end if
        elif (nhbinc  ==  -3) :
c         Possible Lambda helix
          ihfound=1
          nstep=1
          ir=ir+1
          while (ihfound  ==  1)
            nhbinc1=ihbneig(ir)-(ir)
            if (nhbinc1  !=  nhbinc) :
              ihfound=0
            ## end if
c           write (77,6622) ir,nhbinc1,nhbinc2,nhbinc3,ihfound
c6622       format(' ir=',i4,' dihn1,2,3=',3i3,' ifhound=',i2)
            if (ihfound  ==  1) :
              nstep=nstep+1
              irprev=ir
              ir=ir+1
            ## end if
          
          if (ir-irf  >  1) :
            nss=nss+1
            ifss(nss)=irf
            ilss(nss)=ir+nhbinc-1
            itypss(nss)=9
c           write (77,*) 'L Helix nss,ifss,ilss,itypss=',nss,ifss(nss),
          ## end if
        elif (ihbneig(ir)  >  0) :
c         Possible sheet
          isfound=1
          isfirst=0
          islast=0
          npar=0
          napar=0
          ineigmax=0
          ineigmin=nres
          ifs=-1
          while (isfound  ==  1)
            jr=ihbneig(ir)
c           write (77,*) 'Sheet? ir,jr,ihbneig(jr)=',ir,jr,ihbneig(jr)
            if (jr  >  0) :
              if (ihbneig(jr)  ==  ir) :
c               HB(i,j) = HB(j,i) => Antiparallel
                iantiparal(ir)=1
                ing=jr
                ifs=ir-1
                ils=ir+1
c               write (77,*) 'A1 ir=',ir
              ## end if
            ## end if
            if (jr  >  2) :
              if (ihbneig(jr-2)  >  0) :
                if (ihbneig(jr-2)  ==  ir) :
c                 HB(j-1,i)=HB(i,j+1) => Parallel
                  iparal(ir)=1
                  ing=jr
                  ifs=ir
                  ils=ir+2
c                 write (77,*) 'P2 ir=',ir
                ## end if
              ## end if
            ## end if
            if (ir  >  1  and  ir  <  nres) :
              jrm=ihbneig(ir-1)
              if (jrm  >  0) :
                if (ihbneig(jrm)  ==  ir+1) :
c                 HB(i-1,j)=HB(j,i+1) => Parallel
                  iparal(ir)=1
                  ing=jrm
                  ifs=ir
                  ils=ir+2
c                 write (77,*) 'P1 ir=',ir
                ## end if
              ## end if
            else:
              jrm=0
            ## end if
            if (jrm  >  2) :
              if (ihbneig(jrm-2)  >  0) :
                if (ihbneig(jrm-2)  ==  ir+1) :
c                 HB(i-1,j+1)=HB(j-1,i+1) => Antiparallel
                  iantiparal(ir)=1
                  ing=jrm
                  ifs=ir
                  ils=ir+2
c                 write (77,*) 'A2 ir=',ir
                ## end if
              ## end if
            ## end if
            if (iparal(ir)+iantiparal(ir)  >  0) :
              if (isfirst  ==  0) isfirst=ifs
              if (ineigmin  >  ing) ineigmin=ing
              if (ineigmax  <  ing) ineigmax=ing
            elif (ir  ==  1) :
              isfound=0
            elif (ir  >  1) :
              if (iparal(ir-1)+iantiparal(ir-1)  ==  0) :
                isfound=0
                islast=ils
              ## end if
            ## end if
            npar=npar+iparal(ir)
            napar=napar+iantiparal(ir)
            ir=ir+1
c           write (77,*)' isfound,ir=',isfound,ir
          
c         write (77,*)'islast,isfirst,npar,napar,ineigmin,max=',
c    -       islast,isfirst,npar,napar,ineigmin,ineigmax
          if (islast-isfirst  >  2  and  isfirst  >  0) :
            nss=nss+1
            ifss(nss)=isfirst
            ilss(nss)=islast
            if (npar*napar  >  0) :
              itypss(nss)=3
            elif (npar  >  0)  :
              itypss(nss)=1
              if (ineigmax-ineigmin  >  (islast-isfirst)*2)
     -          itypss(nss)=4
            elif (napar  >  0)  :
              itypss(nss)=2
              if (ineigmax-ineigmin  >  (islast-isfirst)*2)
     -          itypss(nss)=5
            ## end if
c           write (77,*) 'Sheet nss,ifss,ilss,itypss=',nss,ifss(nss),
c    -        ilss(nss),itypss(nss)
          ## end if
        ## end if
        if (nss  ==  maxss) :
          write (6,1002) maxss
          if (iwdssp  >  0) write (iwdssp,1002) maxss
          ir=nres
        elif (nss  ==  nss0) :
c         Neither helix nor sheet - just skip over
c         while (ir  <  nres  and  ihbneig(ir)  !=  0)
            ir=ir+1
c         
        else:
          if (nss  >  1) :
            if (ilss(nss-1)  >=  ifss(nss)) ifss(nss)=ilss(nss-1)+1
          else:
            if (ifss(1)  <  1) ifss(1)=1
          ## end if
          ir=ir+1
        ## end if
      
      for iss in range(0, nss):
c        write (77,1001) iss,ifss(iss),ilss(iss),itypss(iss)
c1001    format(' SS',i4,' Start at',i4,' # end at',i4,' Type=',i2)
        for ir in range(ifss(iss), ilss(iss)):
          dssplab(ir)=typc(itypss(iss))
          idistdssp(itypss(iss),ir)=idistdssp(itypss(iss),ir)+1
        
      
      if (iwdssp  >  0) :
        if (nconfig  <=  1) :
          if (nss  >  0) :
            write (iwdssp,2005) (i,ires0+ifss(i),ires0+ilss(i),
     -        ssname(itypss(i))(1:lssname(itypss(i))),i=1,nss)
            write (6,2005) (i,ires0+ifss(i),ires0+ilss(i),
     -        ssname(itypss(i))(1:lssname(itypss(i))),i=1,nss)
          else:
            write (iwdssp,2006)
            write (6,2006)
          ## end if
        ## end if
        if (iwhead  ==  1)
     -     write (iwdssp,2004) (typc(i),ssname(i)(1:lssname(i)),i=1,9)
        iresf=1
        while (iresf  <=  nres)
          iresl=min0(nres,iresf+49)
          for ic in range(0, 50):
            charl1(ic)=' '
            charl2(ic)=' '
          
          for ir in range(max0(3, iresf), min0(iresl, nres-2)):
            if (ixa(ir-2)*ixa(ir)*ixa(ir+2)  !=  0) :
              call angles(dist2(c(1,ixa(ir-2)),c(1,ixa(ir))),
     -          dist2(c(1,ixa(ir+2)),c(1,ixa(ir))),
     -          dist2(c(1,ixa(ir-2)),c(1,ixa(ir+2))),ca1,ca2,cb# end)
              cosa=dble(cb# end)
              b# end=180.0-dacoscheck(cosa,ccc,1,6,'DSSP')*radtodeg
c             write (77,*) ir,' b# end=',b# end
              if (b# end  >  70.0) charl1(ir-iresf+1)='S'
            else:
              charl1(ir-iresf+1)='?'
            ## end if
          
          for ir in range(max0(2, iresf), min0(iresl, nres-2)):
            if (ixa(ir-1)*ixa(ir)*ixa(ir+1)*ixa(ir+2)  !=  0) :
              tors=dihangl(c,ixa(ir-1),ixa(ir),ixa(ir+1),ixa(ir+2),0,n)
c             write (77,*) ir,' tors=',tors*radtodeg
              if (tors  >=  0.0) :
                charl2(ir-iresf+1)='+'
              else:
                charl2(ir-iresf+1)='-'
              ## end if
            else:
              charl2(ir-iresf+1)='?'
            ## end if
          
          write (iwdssp,2000) iresf,iresl,(mod(i,10),i=iresf,iresl)
          write (iwdssp,2001) (mod(ihbneig(i)/1000,10),i=iresf,iresl)
          write (iwdssp,2001) (mod(ihbneig(i)/100,10),i=iresf,iresl)
          write (iwdssp,2001) (mod(ihbneig(i)/10,10),i=iresf,iresl)
          write (iwdssp,2001) (mod(ihbneig(i),10),i=iresf,iresl)
          write (iwdssp,2002) (charl1(i-iresf+1),i=iresf,iresl)
          write (iwdssp,2002) (charl2(i-iresf+1),i=iresf,iresl)
          write (iwdssp,2002) (dssplab(i),i=iresf,iresl)
          iresf=iresl+1
        
      ## end if
      if (iwrose  >  0) :
c       Tentative alternative for turn detection (see Rose's 1977 paper)
        write (iwrose,*) 'Data for turn detection with GW Rose method'
        call normplane(c(1,ixa(1)),c(1,ixa(3)),c(1,ixa(5)),rnprev)
        call normplane(c(1,ixa(2)),c(1,ixa(3)),c(1,ixa(4)),rn1prev)
        for ir in range(4, nres-2):
          call radcirc(c(1,ixa(ir-2)),c(1,ixa(ir)),c(1,ixa(ir+2)),r)
          call normplane(c(1,ixa(ir-2)),c(1,ixa(ir)),c(1,ixa(ir+2)),rn)
          rnn=scprod(rn,rnprev)
          call normplane(c(1,ixa(ir-1)),c(1,ixa(ir)),c(1,ixa(ir+1)),rn1)
          rnn1=scprod(rn1,rn1prev)
          call angles(dist2(c(1,ixa(ir-2)),c(1,ixa(ir))),
     -      dist2(c(1,ixa(ir+2)),c(1,ixa(ir))),
     -      dist2(c(1,ixa(ir-2)),c(1,ixa(ir+2))),ca1,ca2,cb# end)
          cosa=dble(cb# end)
          b# end=180.0-dacoscheck(cosa,ccc,1,6,'DSSP')*radtodeg
          charl1(1)=' '
          if (b# end  >  70.0) charl1(1)='S'
          write (iwrose,2003) ir,charl1(1),dssplab(ir),r,rnn,rnn1
          call trnsfr(rnprev,rn,3)
          call trnsfr(rn1prev,rn1,3)
        
      ## end if
      return
999   write (6,1000) ia,line
      if (iwdssp  >  0) write (iwdssp,1000) ia,line
      return
1000  format(' ERROR: invalid residue number for atom ',i6,':',/,a)
1002  format(' ERROR: maximum number of secondary structure elements (',
     -  i3,') has been reached',/,' - redimension the def ')
1157  format(' eHB update ',i4,' to ',i4,' eold=',f6.2,' old partner=',
     -  i4,' enew=',f6.2)
2000  format[,i6,'-',i5,': ',50i1)
2001  format(14x,50i1)
2002  format(14x,50a1)
2003  format(i5,' DSSP labels:',a1,1x,a1,' rc=',f8.3,
     -  ' rn.rnprv(2)=',f8.4,' rn.rnprv(1)=',f8.4)
2004  format[' line 1: residue number (mod 10)',/,
     -  ' lines 2-5: the digits of the residue number to which the ',/,
     -  12x,'residue number of line 1 is H-bonded (if any)',/,
     -  ' line 6: S for residues with b# end angle ',
     -  '(CA(ir-2)-CA(ir)-CA(ir+2)) > 70 deg',/,
     -  ' line 7: + or -, the sign of the ',
     -  '(CA(ir-1)-CA(ir)-CA(ir+1)-CA(ir+2)) angle',/,
     -  ' line 8: secondary structure element type:',/,
     -  9(9x,a1,': ',a))
2005  format[,(' SS#',i3,' Residue index range: [',i5,',',i5,'] Type:',
     -  a))
2006  format(' No secondary structure element was found')
      # end
